function [ini_tmpl,opt_tmpl,perf] = optim_template(TS,beh,numfr,varargin)
% OPTIM_TEMPLATE      Optimize template based on brain-behavior association
%
%   [ini_tmpl,opt_tmpl,perf] = optim_template(TS,beh,numfr,inds_cv,ini_tmpl,H,hfrac,T0)
%
%   This function takes multiple time series (TS as time x node x subject)
%   and uses simulated annealing to produce a nodal template for filtering
%   time points that improve associations between TS and a behavioral
%   measure (beh). Connectivity components selected and optimized based on
%   frame subsets of size numfr. Each optional input should be called by
%   the parameter name followed by the value or vector input (e.g.,'H',5000
%   or 'inds_cv',[1:2:S]' where S is total number of subjects).
%
%   Inputs:     TS       = time series (time x node x subject x run)
%               beh      = regressed behavior for optimization (subject x 1)
%               numfr    = number of frames filtered by template for optimization
%               inds_cv  = (optional) indices of training subjects
%                                     = (default) order of subjects in TS
%               ini_tmpl = (optional) initial template prior to optimization (time x 1)
%                                     = (default) random filter
%               H        = (optional) number of iterations
%                                     = (default) H = 10000
%               hfrac    = (optional) steepness of temperature gradient (1-hfrac/H)
%                                     = (default) hfrac = 100
%               T0       = (optional) sets initial temperature (and scales the energy term)
%                                     = (default) T0 = 1
%   Outputs:    ini_tmpl = initial template (time x 1)
%               opt_tmpl = optimized template (time x 1)
%               perf     = performance of cost function for current
%                           template and best template at each iteration H

% internal settings
lts = size(TS,1); N = size(TS,2); K = (N^2-N)/2;
S = size(TS,3); F = size(TS,4);
inds = find(triu(ones(N),1));

% TS converted to bipartitions
CIS = false(lts,N,S,F);
for s=1:S
    for f=1:F
        temp = zscore(squeeze(TS(:,:,s,f)));
        CIS(:,:,s,f) = temp>0;
    end
end

% optional parameters
options  = {'inds_cv','ini_tmpl','H','hfrac','T0'};
param = cell(1,length(options));
if nargin > 3
    for i=1:length(options)
        flag = cellfun(@(x) ischar(x) && strcmp(x,options{i}),varargin);
        if any(flag)==true
            temp = find(flag==true);
            param{i} = varargin{temp+1};
        end
    end
end

% default parameters
defaults = {(1:S)',rand(N,1)>0.5,10000,100,1};
for j=1:length(options)
    if ~any(param{j})
        param{j} = defaults{j};
    end
    cmd = sprintf('%s = param{j};',options{j});
    eval(cmd);
end

disp(param)

% internal parameters
S2 = size(inds_cv,1);
Texp = 1-hfrac/H;
Tall = T0.*Texp.^[1:1:H];
initial_cost = 1000;

cost = @(BB) mean(abs(BB));

% select time series of training subjects
TS_train = squeeze(CIS(:,:,inds_cv,:));

% select behavior of training subjects
behav_train = beh(inds_cv);

% initialize optimization
new_tmpl = ini_tmpl;
opt_tmpl = ini_tmpl;
lowcost = initial_cost;
mincost = initial_cost;
perf = zeros(H,2);

% loop through each template iteration
h = 0;
while h<H

    % iteration
    h = h+1;

    % current temperature
    T = T0*Texp^h;

    % select nodes to flip for temporary template
    tind = randperm(N,ceil(abs(randn)));
    if h==1
        tind = []; % preserve initial template
    end
    tmp_tmpl = new_tmpl; tmp_tmpl(tind) = ~new_tmpl(tind);

    % compute MI time series for this temporary template, select top 'numfr'
    % frames, then compute connectivity component
    brain_train = zeros(K,S2);
    % loop over subjects and runs
    parfor s=1:S2
        AGtvec = zeros(K,F);
        for f=1:F
            % compute MI between each frame and the current template
            MI = nmi2_opt2(squeeze(TS_train(:,:,s,f))',tmp_tmpl);
            [~,fr] = maxk(MI,numfr);
            % select frame subset
            cisfr = double(squeeze(TS_train(fr,:,s,f))+1)';
            % compute agreement matrix
            anull = sum(sum([sum(cisfr==1)./N; sum(cisfr==2)./N].*...
                [(sum(cisfr==1)-1)./(N-1); (sum(cisfr==2)-1)./(N-1)],2));
            AG = (agreement2(cisfr) - anull)./numfr;
            AGtvec(:,f) = AG(inds);
        end
        brain_train(:,s) = mean(AGtvec,2);
    end

    % correlate brain and behavior
    [BB, BBp] = corr(brain_train',behav_train,'type','s','rows','complete');

    % cost metric
    costnew = -cost(BB);

    % save current cost
    perf(h,1) = costnew;

    % annealing
    randcon = rand < exp(-(costnew-lowcost)/T);
    if (costnew < lowcost) || randcon
        new_tmpl = tmp_tmpl;
        lowcost = costnew;
        % is this a new absolute best?
        if (lowcost < mincost)
            opt_tmpl = new_tmpl;
            mincost = lowcost;
        end
    end

    % save current best
    perf(h,2) = mincost;

    dispname = ['iteration: ' num2str(h)];
    disp(dispname)
end
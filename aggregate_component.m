function [Xc,moments] = aggregate_component(TS,frames,type,fig,cmap)
% AGGREGATE_COMPONENT     Create connectivity component from frame subset
%
%   [Xc,moments] = aggregate_component(TS,frames,type,fig,cmap)
%
%   This function creates connectivity components for subselections of
%   frames based on input of a frame mask and designation of component type
%   (edge time series (ets) or bipartitions).
%
%   Inputs:     TS      = time series (time x node)
%               frames  = (optional) logical vector(s) of frames selected for bins
%                                    = (default) full time series
%               type    = (optional) method of aggregating frames into component
%                                 1  = (default) averaged from ets
%                                 2  =           aggrement across bipartitions
%               fig     = (optional) create figure of Xc
%                              false = (default) no figure
%                              true  =           figure for each bin in frames
%               cmap    = (optional) colormap (gradient by RGB)
%                                    = (default) parula
%
%   Outputs:    Xc      =  node by node connectivity component matrix
%               moments =  ets or bipartitions used to create Xc

% internal settings
lts = size(TS,1); N = size(TS,2);
inds = find(triu(ones(N),1));

% default parameters
if ~exist('frames','var') || isempty(frames)
    frames = true(lts,1);
end

if ~exist('type','var') || isempty(type)
    type = 1;
end

if ~exist('cmap','var') || isempty(cmap)
    cmap = parula;
end

if ~exist('fig','var') || isempty(fig)
    fig = false;
end

% prepare TS based on type
if type == 1
    ts = double(fcn_edgets(TS));
elseif type == 2
    zTS = zscore(TS);
    ts = double(zTS>0)+1;
end

numfr = sum(frames,1);      % number of frames
numbins = size(frames,2);   % number of bins

% aggregating frames from each bin into connectivity components
for i=1:numbins
    temp = squeeze(frames(:,i));
    bin = squeeze(ts(temp,:))';

    if type==1
        component = mean(bin,2);
    elseif type==2
        anull = sum(sum([sum(bin==1)./N; sum(bin==2)./N].*...
            [(sum(bin==1)-1)./(N-1); (sum(bin==2)-1)./(N-1)],2));
        AG = (agreement2(bin) - anull)./numfr(i);
        component = AG(inds);
    end
    Xc(:,i) = component;
    moments{i} = bin;

    % figure of connectivity component
    if fig == true
        Xcpl = zeros(N);
        Xcpl(inds) = component';
        Xcpl = Xcpl+Xcpl';
        
        figure; imagesc(Xcpl); axis square; colormap(cmap);
        name = ['bin ' num2str(i)]; title(name)
    end
end
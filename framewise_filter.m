function [frames] = framewise_filter(filter,varargin)
% FRAMEWISE_FILTER      Provides frame selection for filter
%
%   [frames] = framewise_filter(filter,numbins,selfr,order,strat)
%
%   This function takes a time series and binning strategy as input and
%   outputs mask(s) for frame selection of each bin. Each optional input
%   should be called by the parameter name followed by the value or option
%   (e.g., 'numbins',12 or 'strat','random').
%
%   Inputs:     filter  = time series of framewise property used for filter                   
%               numbins = (optional) select number of bins
%                                    = (default) 10
%               selfr   = (optional) select vector of bins for output
%                                    = (default) all bins (e.g., 1:numbins)
%               order   = (optional) direction bins are ordered by filter
%                          'descend' = (default) highest to lowest bins
%                          'ascend'  =           lowest to highest bins
%               strat   = (optional) strategy to handle extra frames
%                           false    = (default) discarded from end
%                          'top'     =           distributed over top bins
%                          'bottom'  =           distributed over bottom bins
%                          'random'  =           distributed over random bins
%                          'addfin'  =           larger final bin
%                          'remfin'  =           smaller final bin
%
%   Outputs:    frames  = logical vector(s) of frames selected for bins

% internal settings
options  = {'numbins','selfr','order','strat'};
defaults = {10,1:10,'descend',false};
lts = size(filter,1);

% optional parameters
param = cell(1,length(options));
if nargin > 1
    for i=1:length(options)
        flag = cellfun(@(x) ischar(x) && strcmp(x,options{i}),varargin);
        if any(flag)==true
            temp = find(flag==true);
            param{i} = varargin{temp+1};
        end
    end
end

% default parameters
for j=1:length(options)
    if ~any(param{j})
        if ~any(param{2}) & param{1} ~= defaults{1}
            param{2} = 1:param{1};
        else
            param{j} = defaults{j};
        end
    end
    cmd = sprintf('%s = param{j};',options{j});
    eval(cmd);
end

disp(param)

% assign bin size and frame selection
if strcmp(strat,'remfin')
    numfr = ceil(lts/numbins);
else
    numfr = floor(lts/numbins);
end

fr = repmat(numfr,numbins,1);   % number of frames per bin
x = abs(lts - (numfr*numbins)); % number of extra frames after binning

% reassign extra frames based on strategy
if strat ~= false
    switch strat
        case 'top'
            binsel = 1:x;
        case 'bottom'
            binsel = (param{1}-x)+1:param{1};
        case 'random'
            binsel = randperm(param{1},x);
        case {'addfin','remfin'}
            binsel = param{1};
    end

    distr = length(binsel);
    if distr==1
        addon = x;
    else
        addon = 1;
    end

    if strcmp(strat,'remfin')
        addon = gnegate(addon);
    end

    fr(binsel) = fr(binsel)+addon;
end

% sort frames into each bin
frames = zeros(lts,numbins);
[~,idx] = sort(filter,order);
for bin=1:numbins
    prev = sum(fr(1:(bin-1)));
    inds = idx(prev+1:prev+fr(bin));
    frames(inds,bin) = 1;
end

% create mask(s) for frame selection
frames = logical(squeeze(frames(:,selfr)));
disp(sum(frames))
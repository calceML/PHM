clc, clear all, close all;

data = csvread('.\mvdata.csv');
labeled = '1';
if labeled==1
    nanpos = randi([1 size(data,1)], floor(size(data,1)*0.1), size(data,2)-1);
else
    nanpos = randi([1 size(data,1)], floor(size(data,1)*0.1), size(data,2));
end

for m=1 : size(nanpos,2)
    data(nanpos(:,m),m) = NaN;
end

%% Missing value management with no moving window
% dataimpute = calcemvmangement(data, labeled, 'MissingValManMovMethod', 'movmedian', 'MovWindows', '20',...
%     'PlotXDim', '1', 'PlotYDim', '2');

%% Missing value management with moving window
dataimpute = calcemvmangement(data, labeled, 'MissingValManMovMethod', 'movmedian', 'MissingValManMethod', 'linear', 'MovWindows', '20',...
    'PlotXDim', '1', 'PlotYDim', '2');

%% Missing Value Management
function [dataimpute] = calcemvmangement(data, labeled, varargin)

p = inputParser;
p.CaseSensitive = false;    % Names are not sensitive to case: 'a' matches 'A'

defaultPlotXDim = '1';
defaultPlotYDim = '2';
defaultMissingValManMethod = 'nearest'; % methods ={'linear', 'nearest'}
defaultMissingValManMovMethod = 'none'; % movmethods = {'movmean', 'movmedian'}
defaultMovWindows = '5';

addRequired(p, 'data', @ismatrix);
addRequired(p, 'labeled', @ischar);
addParameter(p, 'plotxdim', defaultPlotXDim, @ischar);
addParameter(p, 'plotydim', defaultPlotYDim, @ischar);
addParameter(p, 'missingvalmanmethod', defaultMissingValManMethod, @ischar);
addParameter(p, 'missingvalmanmovmethod', defaultMissingValManMovMethod, @ischar);
addParameter(p, 'movwindows', defaultMovWindows, @ischar);
parse(p, data, labeled, varargin{:});

labeled = str2num(p.Results.labeled);
if labeled==1
    data = p.Results.data(:,1:end-1);
    label = p.Results.data(:,end);
else
    data = p.Results.data;
end

if strcmp(p.Results.missingvalmanmovmethod, 'none')
    ndims = size(data,2);
    dataimpute = [];
    nanidx = [];
    for m=1 : ndims
        tmp = data(:,m);
        nanidx(:,m) = isnan(tmp);
        if strcmp(p.Results.missingvalmanmethod, 'linear')
            dataimpute(:,m) = fillmissing(tmp, 'linear');
        elseif strcmp(p.Results.missingvalmanmethod, 'nearest')
            dataimpute(:,m) = fillmissing(tmp, 'linear');
        end
    end
    
    if ndims==1
        calcemvmplot(dataimpute, nanidx, [str2num(p.Results.plotxdim)]);
    else
        calcemvmplot(dataimpute, nanidx, [str2num(p.Results.plotxdim) str2num(p.Results.plotydim)]);
    end
    
    if labeled==1
        dataimpute = [dataimpute label];
    end
    
elseif strcmp(p.Results.missingvalmanmovmethod, 'movmean')
    ndims = size(data,2);
    dataimpute = [];
    nanidx = [];
    for m=1 : ndims
        tmp = data(:,m);
        nanidx(:,m) = isnan(tmp);
        dataimpute(:,m) = fillmissing(tmp, 'movmean', str2num(p.Results.movwindows));
    end
    
    if ndims==1
        calcemvmplot(dataimpute, nanidx, [str2num(p.Results.plotxdim)]);
    else
        calcemvmplot(dataimpute, nanidx, [str2num(p.Results.plotxdim) str2num(p.Results.plotydim)]);
    end
    
    if labeled==1
        dataimpute = [dataimpute label];
    end
elseif strcmp(p.Results.missingvalmanmovmethod, 'movmedian')
    ndims = size(data,2);
    dataimpute = [];
    nanidx = [];
    for m=1 : ndims
        tmp = data(:,m);
        nanidx(:,m) = isnan(tmp);
        dataimpute(:,m) = fillmissing(tmp, 'movmedian', str2num(p.Results.movwindows));
    end
    
    if ndims==1
        calcemvmplot(dataimpute, nanidx, [str2num(p.Results.plotxdim)]);
    else
        calcemvmplot(dataimpute, nanidx, [str2num(p.Results.plotxdim) str2num(p.Results.plotydim)]);
    end
    
    if labeled==1
        dataimpute = [dataimpute label];
    end
end
    
end

function calcemvmplot(dataimpute, nanidx, plotdim)

nanidx = logical(nanidx);
ndim = length(plotdim);

if ndim==1
    tmpaxis = [1:size(dataimpute,1)]';
    numnan = sum(nanidx);
    
    figure, grid on;
    plot(tmpaxis, dataimpute, 'or', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
    hold on, plot(tmpaxis(nanidx), dataimpute(nanidx,plotdim), 'sb', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
    xlabel('# of samples'), ylabel(['Dim: ' num2str(plotdim)]);
    legend(['Raw Data', 'Imputed Data']);
    title(['Missting Value Management Results' '(# of NaNs = ' num2str(numnan) ' )']);
    
elseif ndim>=2
    
    tmpaxis = [1:size(dataimpute,1)]';
    numnan = sum(nanidx(:,plotdim));
    
    figure, grid on;
        
    subplot(2,1,1), plot(tmpaxis, dataimpute(:,plotdim(2)),...
        'ro', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
    hold on, plot(tmpaxis(nanidx(:,plotdim(2))), dataimpute(nanidx(:,plotdim(2)),plotdim(2)),...
        'bo', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
    xlabel('# of samples'), ylabel(['Dim: ' num2str(plotdim(2))]);
    legend('Raw Data', 'Imputed Data');
    title(['Missting Value Management Results' ' (# of Initial NaNs = ' num2str(numnan(2)) ' )']);
    
    subplot(2,1,2), plot(tmpaxis, dataimpute(:,plotdim(1)),...
        'ro', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
    hold on, plot(tmpaxis(nanidx(:,plotdim(1))), dataimpute(nanidx(:,plotdim(1)),plotdim(1)),...
        'bo', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
    xlabel('# of samples'), ylabel(['Dim: ' num2str(plotdim(1))]);
    legend('Raw Data', 'Imputed Data');
    title(['Missting Value Management Results' ' (# of Initial NaNs = ' num2str(numnan(1)) ' )']);
    
    
end
    
end
function calceocsvmtraining(data, varargin)

% [ocsvmmodel] = calceocsvmtraining(data, 'IsLabel', 'true', 'EnablePCA', 'false', 'KernelMethod', 'linear', 'OurlierFraction', '0.05', 'PlotXDim', '1', 'PlotYDim', '2');

p = inputParser;
p.CaseSensitive = false;    % Names are not sensitive to case: 'a' matches 'A'

defaultIsLabel = 'true';
defaultEnablePCA = 'false';
defaultKernelMethod = 'linear';
expectedKernelMethods = {'linear', 'polynomial', 'gaussian'};
defaultPolynomialOrder = '2';
defaultGaussianSigma = '1';
defaultOutlierFraction = '0.2';
defaultPlotXDim = 'Dim1';
defaultPlotYDim = 'Dim2';

addRequired(p, 'data', @ismatrix);
addParameter(p, 'islabel', defaultIsLabel, @ischar);
addParameter(p, 'enablepca', defaultEnablePCA, @ischar);
addParameter(p, 'kernelmethod', defaultKernelMethod,...
    @(x) any(validatestring(x, expectedKernelMethods)));
addParameter(p, 'polynomialorder', defaultPolynomialOrder, @ischar);
addParameter(p, 'gaussiansigma', defaultGaussianSigma, @ischar);
addParameter(p, 'outlierfraction', defaultOutlierFraction, @ischar);
addParameter(p, 'plotxdim', defaultPlotXDim, @ischar);
addParameter(p, 'plotydim', defaultPlotYDim, @ischar);

parse(p, data, varargin{:});

%% To configure data & label

if strcmp(p.Results.islabel, 'true')
    x = data(:,1:end-1);
    y = data(:,end);
    if strcmp(p.Results.enablepca, 'true')
        % need to load pc
        x = x*pc;
    end
    
    if size(x,2) < 2
        plotxdim = 1;
        plotydim = [];
    else
        plotxdim = str2num(p.Results.plotxdim(isstrprop(p.Results.plotxdim, 'digit')));
        plotydim = str2num(p.Results.plotydim(isstrprop(p.Results.plotydim, 'digit')));
    end
    
elseif strcmp(p.Results.islabel, 'false')
    x = data;
    y = zeros(size(data,1),1);
    if strcmp(p.Results.enablepca, 'true')
        % need to load pc
        x = x*pc;
    end
    
    if size(x,2) < 2
        plotxdim = 1;
        plotydim = [];
    else
        plotxdim = str2num(p.Results.plotxdim(isstrprop(p.Results.plotxdim, 'digit')));
        plotydim = str2num(p.Results.plotydim(isstrprop(p.Results.plotydim, 'digit')));
    end
    
end

%% To train an ocsvm model

rng(1);

if strcmp(p.Results.kernelmethod, 'linear')
    ocsvmmodel = fitcsvm(x(:, [plotxdim, plotydim]), y, 'Standardize', true, 'KernelFunction', 'linear', 'BoxConstraint', 1, 'KernelScale', 'auto',...
        'OutlierFraction', str2num(p.Results.outlierfraction));
elseif strcmp(p.Results.kernelmethod, 'polynomial')
    ocsvmmodel = fitcsvm(x(:, [plotxdim, plotydim]), y, 'Standardize', true, 'KernelFunction', 'polynomial', 'PolynomialOrder', str2num(p.Results.polynomialorder),...
        'BoxConstraint', 1, 'KernelScale', 'auto', 'OutlierFraction', str2num(p.Results.outlierfraction));
elseif strcmp(p.Results.kernelmethod, 'gaussian')
    ocsvmmodel = fitcsvm(x(:, [plotxdim, plotydim]), y, 'Standardize', true, 'KernelFunction', 'gaussian', 'BoxConstraint', 1,...
        'KernelScale', str2num(p.Results.gaussiansigma)', 'OutlierFraction', str2num(p.Results.outlierfraction));
end

saveFolder = 'C:\CALCE-PHM\OCSVMAD';
if ~exist(saveFolder, 'dir');
    mkdir(saveFolder);
end
save([saveFolder '\ocsvmmodel.mat'], 'ocsvmmodel');

%% To plot the result of ocsvm training..

calceocsvmtrainingplot(x, plotxdim, plotydim);




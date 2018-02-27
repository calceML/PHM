function calceocsvmtest(data, varargin)

p = inputParser;
p.CaseSensitive = false;    % Names are not sensitive to case: 'a' matches 'A'

defaultIsLabel = 'true';
defaultEnablePCA = 'false';
defaultPlotXDim = 'Dim1';
defaultPlotYDim = 'Dim2';

addRequired(p, 'data', @ismatrix);
addParameter(p, 'islabel', defaultIsLabel, @ischar);
addParameter(p, 'enablepca', defaultEnablePCA, @ischar);
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

%% To plot anomaly detection result..

calceocsvmtestplot(x, plotxdim, plotydim);

end



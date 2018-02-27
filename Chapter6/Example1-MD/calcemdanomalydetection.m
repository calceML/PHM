function calcemdanomalydetection(refdata, data, varargin)

p = inputParser;
p.CaseSensitive = false;    % Names are not sensitive to case: 'a' matches 'A'

defaultThreshold = '';

addRequired(p, 'refdata', @ismatrix);
addRequired(p, 'data', @ismatrix);
addParameter(p, 'threshold', defaultThreshold, @ischar);
parse(p, refdata, data, varargin{:});

mu = mean(refdata);
SIGMA = cov(refdata);
mdrefdata = calcemd(refdata, mu, SIGMA);

threshold = str2num(p.Results.threshold);
if isempty(threshold)
    threshold = mean(mdrefdata) + 3*std(mdrefdata);
end

mddata = calcemd(data, mu, SIGMA);

%% to plot md values
calcemdplot(mdrefdata, mddata, threshold);

end



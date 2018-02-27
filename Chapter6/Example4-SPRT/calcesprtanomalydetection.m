function [sprtindex, alarm, A, B] = calcesprtanomalydetection(refdata, testdata, varargin)

% This function applies the Sequential Probability Ratio Test (SPRT) for anomaly detection.
%
% Input:
%   refdata: a column vector involving reference data (healthy data)
%   data: a column vector involving samples that are tested for anomaly detection
%   alpha: false-alarm probability, ranging from 0 to 1
%   beta: missed-alarm probability, ranging from 0 to 1
%
% Output:
%   sprtindex: SPRT index
%   alarm: if alarm is 1, the data point is detected as an anomaly
%   A, B are lower and higher detection thresholds, respectively.

p = inputParser;
p.CaseSensitive = false;    % Names are not sensitive to case: 'a' matches 'A'

defaultIsLabel = 'true';
defaultSPRTMethod = 'mean';
expectedSPRTMethods = {'mean', 'variance'};
% a pre-assigned scalar used as the system disturbance magnitude for the
% mean test
defaultSystDistMagMean = ''; % 4*std of the reference data
% a pre-assigned scalar used as the system disturbance magnitude for the
% variance test
defaultSystDistMagVariance = '2';
defaultAlpha = '0.025';
defaultBeta = '0.025';

addRequired(p, 'refdata', @ismatrix);
addRequired(p, 'testdata', @ismatrix);
addParameter(p, 'islabel', defaultIsLabel, @ischar);
addParameter(p, 'sprtmethod', defaultSPRTMethod,...
    @(x) any(validatestring(x, expectedSPRTMethods)));
addParameter(p, 'systdistmagmean', defaultSystDistMagMean, @ischar);
addParameter(p, 'systdistmagvariance', defaultSystDistMagVariance, @ischar);
addParameter(p, 'alpha', defaultAlpha, @ischar);
addParameter(p, 'beta', defaultBeta, @ischar);
parse(p, refdata, testdata, varargin{:});

if strcmp(p.Results.islabel, 'true')
    orefdata = refdata;
    otestdata = testdata;
    refdata = refdata(:,1:end-1);
    testdata = testdata(:, 1:end-1);
elseif strcmp(p.Results.islabel, 'false')
    orefdata = refdata;
    otestdata = testdata;
    refdata = refdata;
    testdata = testdata;
end

alpha = str2num(p.Results.alpha);
beta = str2num(p.Results.beta);
A = log(beta/(1-alpha));
B = log((1-beta)/alpha);

systdistmagmean = str2num(p.Results.systdistmagmean);
systdistmagvariance = str2num(p.Results.systdistmagvariance);
if isempty(systdistmagmean)
    systdistmagmean = 4*std(refdata);
end

[mu, ~] = normfit(refdata);
refdata = refdata - mu;
testdata = testdata - mu;
[~, sigma] = normfit(refdata);

index = 0;
accumulation = 0;
nobservations = length(testdata);
sprtindex = zeros(nobservations, 1);
alarm = zeros(nobservations, 1);
if strcmp(p.Results.sprtmethod, 'mean')
    for m=1 : nobservations
        index = index + ((systdistmagmean/(sigma^2))*(testdata(m) - systdistmagmean/2));
        sprtindex(m) = index;
        if index<A
            alarm(m) = 0;
            index = 0;
        elseif index>B
            alarm(m) = 1;
            index = 0;
        end
    end
elseif strcmp(p.Results.sprtmethod, 'variance')
    for m=1 : nobservations
        accumulation = accumulation + (testdata(m)^2);
        index = ((accumulation/(2*(sigma^2)))*(1-(1/systdistmagvariance)))-((m/2)*log(systdistmagvariance));
        sprtindex(m) = index;
        if index<A
            alarm(m) = 0;
        elseif index>B
            alarm(m) = 1;
            accumulation = 0;
        end
    end
end

calcesprtplot(otestdata, alarm, p.Results.islabel);

end
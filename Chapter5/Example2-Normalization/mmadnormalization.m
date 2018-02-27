function output = mmadnormalization(data)

% The median and median absolute deviation (MMAD) normalization method is a
% robust measure of the variability of a univariate sample of quantitative
% data. Median absolute deviation (MAD) is a measure of statistical
% dispersion and ore resilient to outliers in a data set than the standard
% deviation.

if (nargin ~= 1)
    msg = 'Please use a correct input for the MMAD normalization';
    error(msg);
end

mad = median(abs(data - median(data)));
output = (data - median(data))/mad;

end
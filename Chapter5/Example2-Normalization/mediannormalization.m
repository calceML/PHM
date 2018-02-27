function output = mediannormalization(data)

% The median method normalizes each sample by the median of raw inputs for
% all the inputs in the sample. It is a useful normalization to use when
% there is a need to compute the ratio between two hybridized samples.
% Median is not influenced by the magnitude of extreme deviations. It can
% be more useful when performing the distribution.

% Input:
% data is a a x b matrix, where a is the number of data in the ith variable
% and b is the number of variables.

if (nargin ~=1)
    msg = 'Please use an correct input for the median normalization';
    error(msg);
end

a = size(data,1);
b = size(data,2);
output = zeros(a,b);
for i=1 : b
    output(:,i) = data(:,i)/median(data(:,i));
end

end
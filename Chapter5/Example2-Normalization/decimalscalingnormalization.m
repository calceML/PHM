function output = decimalscalingnormalization(data)

% In the decimal scaling normalization, the decimal point of the values of
% an attribute 'data' is moved to its maximum absolute value. The number of
% decimal points moved depends on the maximum absolute value of the data
% set 'data.'

if (nargin ~= 1)
    msg = 'Please use correct inputs for the decimal scaling normalization';
    error(msg);
end

ndecpnt = ceil(log10(max(abs(data))));
output = data/10^ndecpnt;

end
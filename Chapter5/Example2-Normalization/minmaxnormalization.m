function output = minmaxnormalization(data, minbound, maxbound)

% Min-Max Normalization is a simple technique where the technique can
% specifically fit the data in a pre-defined boundary with a pre-defined
% boundary.
% inputs:
% data should be a 1 x a vector, where a is the number of data for a certain variable.
% [minbound maxbound] is the pre-defined boundary. By default, minbound is
% 0 and maxbound is 1.
% ouput:
% output is the min-max normalized data

if (nargin > 3 | nargin < 1)
    msg = 'Please use an correct input(s) for the min-max normalization';
    error(msg);
elseif (nargin==1)
    minbound = 0;
    maxbound = 1;
end

minvalue = min(data);
maxvalue = max(data);
n = data - minvalue;
d = maxvalue - minvalue;
output = (n/d)*(maxbound - minbound) + minbound;
    
end
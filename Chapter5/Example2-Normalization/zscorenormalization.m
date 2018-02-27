function output = zscorenormalization(data)

% Z-Score Normalization gives the normalized values or range of data from
% the original unstructured data using the concepts like mean and standard
% deviation. That is, the unstructured data can be normalized using z-score
% parameter
% input:
% data is a a x b matrix, where a is the dimensionality of ith variable and
% b is the number of variables.
% output is z-score normalized values.

a = size(data,1);
b = size(data,2);
output = zeros(a,b);
for i=1 : b
    m = mean(data(:,i));
    s = std(data(:,i));
    output(:,i) = (data(:,i) - m)/s;
end

end
    
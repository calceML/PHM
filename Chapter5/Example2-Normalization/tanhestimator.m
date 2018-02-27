function output = tanhestimator(data)

% Tanh estimator is a robust and highly efficient method to normalize time
% series data.

m = mean(data);
s = std(data);

output = 0.5*(tanh((0.01*(data - m))/s)+1);

end
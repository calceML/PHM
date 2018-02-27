clc, clear all, close all;

data = csvread('timeseries.csv');
data = data(:,2);

method = 'tanh';

switch method
    case 'minmax'
        minbound = 0;
        maxbound = 1;
        normalized = minmaxnormalization(data, minbound, maxbound);
    case 'zscore'
        normalized = zscorenormalization(data);
    case 'decimalscaling'
        normalized = decimalscalingnormalization(data);
    case 'median'
        normalized = mediannormalization(data);
    case 'mmad'
        normalized = mmadnormalization(data);
    case 'tanh'
        normalized = tanhestimator(data);
end

figure,
subplot(2,1,1), plot(data, 'b-');
xlabel('# of cycles'), ylabel('Normalized Discharge Capacity [Ah]');
title('Raw Data');
subplot(2,1,2), plot(normalized, 'r-');
xlabel('# of cycles'), ylabel('Normalized Discharge Capacity [Ah]');
title('Normalized Data');

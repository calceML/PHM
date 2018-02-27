clc, clear all, close all;

% To load a time-series data
tsdata = csvread('timeseries.csv');
predictor = tsdata(:,1);    % # of cycle for a lithium-ion battery cycling test
response = tsdata(:,2);     % discharge capacity associated with the number of cycle

figure, plot(predictor, response, 'b-');
xlabel('# of cycles'), ylabel('Normalized Discharge Capacity [Ah]');
title('Discharge Capacity over Cycles');

% To create a linear regression model
w = polyfit(predictor, response, 1);
estimated = w(1)*predictor + w(2);
residual = estimated - response;

% To plot the results
figure,
subplot(2,1,1), plot(predictor, response, 'b-');
hold on, plot(predictor, estimated, 'r-.');
xlabel('# of cycles'), ylabel('Normalized Discharge Capacity [Ah]');
legend('True', 'Estimated');
title('Linear Regression for Estiamting Discharge Capacity over Cycles');

subplot(2,1,2), plot(predictor, residual, 'g-');
xlabel('# of cycles'), ylabel('Normalized Discharge Capacity [Ah]');
title('Residual');



    

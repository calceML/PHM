clc, clear all, close all;

% To load a time-series data
tsdata = csvread('timeseries.csv');
predictor = tsdata(:,1);    % # of cycle for a lithium-ion battery cycling test
response = tsdata(:,2);     % discharge capacity associated with the number of cycle

figure, plot(predictor, response, 'b-');
xlabel('# of cycles'), ylabel('Normalized Discharge Capacity [Ah]');
title('Discharge Capacity over Cycles');

% To create an elastic net regression model
dist = 'normal';    % dist represents a distributional family for the nonsystematic variation in the responses.
alpha = 0.5;  % alpha from 0 to 1. If alpha=1, it is a LASSO regression model. Otherwise, it is an elastic net regression model. If alpha is coloe to 0, it is a ridge regression model.
[w, fitinfo] = lassoglm(predictor, response, dist, 'Alpha', alpha, 'CV', 5);

% To plot the LASSO regression result
lassoPlot(w, fitinfo, 'PlotType', 'CV');




    

clc, clear all, close all;

%% load data

data = csvread('dataset.csv', 2, 0);
nolabeldata = data(:,1:end-1);
pc1 = calcepca(data, 'IsLabel', 'true', 'HealthyDataOnly', 'false', 'AllDataWithLabel', 'true', 'ExplainedVariance', '95', 'FirstNPCs', '', 'PlotXDim', 'PC1', 'PlotYDim', 'PC2');
% pc2 = calcepca(data, 'IsLabel', true, 'HealthyDataOnly', false, 'AllDataWithLabel', true, 'ExplainedVariance', '', 'FirstNPCs', '4', 'PlotXDim', 'PC2', 'PlotYDim', 'PC3');
% pc3 = calcepca(nolabeldata, 'IsLabel', false, 'HealthyDataOnly', false, 'AllDataWithLabel', false, 'ExplainedVariance', '', 'FirstNPCs', '10', 'PlotXDim', 'PC5', 'PlotYDim', 'PC6');



clc, clear all, close all;

%% load data

data = csvread('AlternatorDataWithLabels.csv', 2, 0);
ldc1 = calcelda(data, 'PlotXDim', '', 'PlotYDim', '');


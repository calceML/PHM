clc, clear all, close all;

refdata = csvread('refdata.csv');
testdata = csvread('testdata.csv');

refdata1 = refdata(:,1);
testdata1 = testdata(:,1);

% [sprtindex, alarm, A, B] = calcesprtanomalydetection(refdata, testdata, 'IsLabel', 'true', 'SPRTMethod', 'variance',...
%     'SystDistMagVariance', '1', 'Alpha', '0.01', 'Beta', '0.01');

[sprtindex, alarm, A, B] = calcesprtanomalydetection(refdata, testdata, 'IsLabel', 'true', 'SPRTMethod', 'mean',...
    'SystDistMagMean', num2str(4*std(refdata(:,1))), 'Alpha', '0.01', 'Beta', '0.01');

% [sprtindex, alarm, A, B] = calcesprtanomalydetection(refdata1, testdata1, 'IsLabel', 'false', 'SPRTMethod', 'mean',...
%     'SystDistMagMean', num2str(4*std(refdata(:,1))), 'Alpha', '0.005', 'Beta', '0.001');
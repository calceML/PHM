clc, clear all, close all

%% load data

data = csvread('dataset.csv', 2, 0);
ndimensions = size(data,2);

refdata = data(data(:,ndimensions)==0,1:ndimensions-1);
tdata = data(data(:,ndimensions)==1,1:ndimensions-1);

calcemdanomalydetection(refdata, tdata, 'Threshold', '');


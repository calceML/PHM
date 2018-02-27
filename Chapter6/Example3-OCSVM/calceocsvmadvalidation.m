clc, clear all, close all;

data = csvread('dataset.csv', 2, 0);
ndimensions = size(data,2);

refdata = data(data(:,ndimensions)==0,1:ndimensions-1);
tdata = data(data(:,ndimensions)==1,1:ndimensions-1);

calceocsvmtraining(refdata, 'IsLabel', 'false', 'EnablePCA', 'false', 'KernelMethod', 'gaussian', 'GaussianSigma', '0.5', 'PlotXDim', 'Dim1', 'PlotYDim', 'Dim2');
calceocsvmtest(tdata, 'IsLabel', 'false', 'EnablePCA', 'false', 'PlotXDim', 'Dim1', 'PlotYDim', 'Dim2');

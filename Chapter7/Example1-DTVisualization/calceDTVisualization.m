clc, clear all, close all;

mcdata = csvread('multiclassdata.csv');

meas = mcdata(:,1:end-1);
labels = mcdata(:,end);

DTMdl = fitctree (meas, labels);
view(DTMdl, 'Mode', 'graph')

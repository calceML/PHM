
clc, clear all, close all;

Md.Name = 'SVMR';
Md.SVMR.KernelName = 'Gaussian';
Md.SVMR.KernelParameter = '2^-3';

% SVMR
Md.SVMR.LabelStr = {'Cycles', 'Normalized Discharge Capacity [Ah]'};
Md.SVMR.LegendStr = {'Response (Test)', 'Predicted Response'};

% data: a M x N matrix, where M is # of observations and N is # of
%       predictors.

data = csvread('timeseries.csv');

switch Md.Name
    case 'SVMR'

        % x: predictors
        % y: response
        x = data(:,1:end-1);
        y = data(:,end);
        if strcmp(Md.SVMR.KernelName, 'Linear')
            Md.SVMR.TrainedModel = fitrsvm(x, y, 'Standardize',true,...
                                'KernelFunction', 'linear',...
                                'BoxConstraint', 1, 'KernelScale', 'auto');
        elseif strcmp(Md.SVMR.KernelName, 'Polynomial')
            Md.SVMR.TrainedModel = fitrsvm(x, y, 'Standardize', true,...
                                'KernelFunction', 'polynomial',...
                                'PolynomialOrder', str2num(Md.SVMR.KernelParameter),...
                                'BoxConstraint', 1, 'KernelScale', 'auto');
        elseif strcmp(Md.SVMR.KernelName, 'Gaussian')
            Md.SVMR.TrainedModel = fitrsvm(x, y, 'Standardize', true,...
                                'KernelFunction', 'gaussian',...
                                'BoxConstraint', 1, 'KernelScale', str2num(Md.SVMR.KernelParameter));
        end

        % Assume testdata is a new test data
        testdata = data;
        testx = testdata(:,1:end-1);
        testy = testdata(:,end);
        yhat = predict(Md.SVMR.TrainedModel, testx);

        % To compute root mean squared error (rmse)
        rmse =sqrt(mean((testy - yhat).^2));

        calcePlotPrognosis(testx, testy, yhat, rmse, Md.SVMR.LabelStr, Md.SVMR.LegendStr);       

    otherwise
        disp('Unknown Prognosis Methods!!');
end

function calcePlotPrognosis(testx, testy, yhat, rmse, labelStr, legendStr)
    figure,
    plot(testy, 'b-');
    hold on, plot(yhat, 'r-.');
    xlabel(labelStr{1});
    ylabel(labelStr{2});

    str = {'RMSE = ', num2str(rmse)};
    fxpos = ceil(median(find(testx>0.8*size(testx,1))));
    fypos = testy(ceil(median(find(testy>0.7*max(testy)))));
    text(fxpos, fypos, str, 'Color', 'red', 'FontSize', 14);
    legend(legendStr);
end







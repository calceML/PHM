clc, clear all, close all;

fcn = 'gaussian';

nClass = 2;

switch fcn
    case 'gaussian'
        %generate random multivariate Gaussian data
        A1 = 0.2*randn(2,2);
        A2 = 0.15*rand(2,2);
        nClass1 = 400;
        nClass2 = 350;
        m1 = 3.4; m2 = 4;
        Y = mvnrnd(0.1*randn(1,2)+m1,(A1*A1'),nClass1);
        Y = [Y;mvnrnd(0.1*randn(1,2)+m2,(A2*A2'),nClass2)];
        labels = [ones(nClass1,1); 2*ones(nClass2,1)]; %create label vector
        if nClass > 2
            nClass3 = 375;
            m3 = 2.6;
            A3 = 0.3*rand(2,2);
            Y = [Y;mvnrnd(0.1*rand(1,2) + m3, A3*A3',nClass3)]; 
            labels = [ones(nClass1,1); 2*ones(nClass2,1);...
                        3*ones(nClass3,1)]; %create label vector
        end

        figure(1);
        for i=1:nClass
            plot(Y(labels==i,1),Y(labels==i,2),'.','MarkerSize',12)
            hold on;
        end

        model = naiveBayesFit(Y,labels,fcn);
        X = Y;

        %create different simulated data from same distributions
        %generate new multivariate Gaussian data for testing
        nClass1 = 200;
        nClass2 = 250;
        Y = mvnrnd(0.1*randn(1,2)+m1+0.05,(A1*A1'),nClass1);
        Y = [Y;mvnrnd(0.1*randn(1,2)+m2-0.05,(A2*A2'),nClass2)];
        labels = [ones(nClass1,1); 2*ones(nClass2,1)]; %create label vector
        if nClass > 2
            nClass3 = 375;
            Y = [Y;mvnrnd(0.1*rand(1,2) + m3-0.05, A3*A3',nClass3)]; 
            labels = [ones(nClass1,1); 2*ones(nClass2,1);...
                        3*ones(nClass3,1)]; %create label vector
        end
        [yhat,P] = naiveBayesPredict(Y,model);

        figure(1);
        for i=1:nClass
            plot(Y(labels==i,1),Y(labels==i,2),'x','MarkerSize',6)
            hold on;
        end
        misLabel = labels~=yhat;
        plot(Y(misLabel,1),Y(misLabel,2),'ko','MarkerSize',8)
        if nClass==2
            legend('Class 1','Class 2','Class 1 Test','Class 2 Test','Misclassified')
        elseif nClass==3
            legend('Class 1','Class 2','Class 3',...
                        'Class 1 Test','Class 2 Test','Class 3 Test','Misclassified')
        end
        plotClassBoundary([X;Y],labels,@(X)naiveBayesPredict(X,model),false);
        set(gca,'FontSize',18)

        C = confusionmat(labels, yhat);
        plotConfMat(C);
 
    case 'lognormal'
        %generate random multivariate lognormal data
        A1 = 0.1*randn(2,2);
        A2 = 0.1*rand(2,2);
        nClass1 = 400;
        nClass2 = 350;
        Y = mvnrnd(0.1*randn(1,2)+3,(A1*A1'),nClass1);
        Y = [Y;mvnrnd(0.025*randn(1,2)+4,(A2*A2'),nClass2)];
        labels = [ones(nClass1,1); 2*ones(nClass2,1)]; %create label vector
        if nClass > 2
            nClass3 = 375;
            A3 = 0.2*rand(2,2);
            Y = [Y;mvnrnd(0.1*rand(1,2)+3.5, A3*A3',nClass3)]; 
            labels = [ones(nClass1,1); 2*ones(nClass2,1);...
                        3*ones(nClass3,1)]; %create label vector
        end
        Y = exp(Y); %transform to log-normal
        figure(1);
        for i=1:nClass
            plot(Y(labels==i,1),Y(labels==i,2),'.','MarkerSize',12)
            hold on;
        end

        model = naiveBayesFit(Y,labels,fcn);
        X = Y;

        %create different simulated data from same distributions
        nClass1 = 200;
        nClass2 = 250;
        Y = mvnrnd(0.1*randn(1,2)+3.05,(A1*A1'),nClass1);
        Y = [Y;mvnrnd(0.025*randn(1,2)+3.98,(A2*A2'),nClass2)];
        labels = [ones(nClass1,1); 2*ones(nClass2,1)]; %create label vector
        if nClass > 2
            nClass3 = 375;
            Y = [Y;mvnrnd(0.1*rand(1,2)+3.45, A3*A3',nClass3)]; 
            labels = [ones(nClass1,1); 2*ones(nClass2,1);...
                        3*ones(nClass3,1)]; %create label vector
        end
        Y = exp(Y);
        [yhat,P] = naiveBayesPredict(Y,model);

        figure(1);
        for i=1:nClass
            plot(Y(labels==i,1),Y(labels==i,2),'x','MarkerSize',6)
            hold on;
        end
        misLabel = labels~=yhat;
        plot(Y(misLabel,1),Y(misLabel,2),'ko','MarkerSize',8)
        if nClass==2
            legend('Class 1','Class 2','Class 1 Test','Class 2 Test','Misclassified')
        elseif nClass==3
            legend('Class 1','Class 2','Class 3',...
                        'Class 1 Test','Class 2 Test','Class 3 Test','Misclassified')
        end
        plotClassBoundary([X;Y],labels,@(X)naiveBayesPredict(X,model),false);
        set(gca,'FontSize',18)

        C = confusionmat(labels, yhat);
        plotConfMat(C);
        
    case 'weibull'
        %generate random multivariate Gaussian data
        A1 = 0.1*randn(2,2);
        A2 = 0.1*rand(2,2);
        nClass1 = 400;
        nClass2 = 350;
        Y = mvnrnd(0.1*randn(1,2)+3,(A1*A1'),nClass1);
        Y = [Y;mvnrnd(0.025*randn(1,2)+4,(A2*A2'),nClass2)];
        labels = [ones(nClass1,1); 2*ones(nClass2,1)]; %create label vector
        if nClass > 2
            nClass3 = 375;
            A3 = 0.2*rand(2,2);
            Y = [Y;mvnrnd(0.1*rand(1,2)+3.5, A3*A3',nClass3)]; 
            labels = [ones(nClass1,1); 2*ones(nClass2,1);...
                        3*ones(nClass3,1)]; %create label vector
        end
        Y = exp(Y); %transform to log-normal
        figure(1);
        for i=1:nClass
            plot(Y(labels==i,1),Y(labels==i,2),'.','MarkerSize',12)
            hold on;
        end

        model = naiveBayesFit(Y,labels,fcn);
        X = Y;

        %create different simulated data from same distributions
        nClass1 = 200;
        nClass2 = 250;
        Y = mvnrnd(0.1*randn(1,2)+3.05,(A1*A1'),nClass1);
        Y = [Y;mvnrnd(0.025*randn(1,2)+3.98,(A2*A2'),nClass2)];
        labels = [ones(nClass1,1); 2*ones(nClass2,1)]; %create label vector
        if nClass > 2
            nClass3 = 375;
            Y = [Y;mvnrnd(0.1*rand(1,2)+3.45, A3*A3',nClass3)]; 
            labels = [ones(nClass1,1); 2*ones(nClass2,1);...
                        3*ones(nClass3,1)]; %create label vector
        end
        Y = exp(Y);
        [yhat,P] = naiveBayesPredict(Y,model);

        figure(1);
        for i=1:nClass
            plot(Y(labels==i,1),Y(labels==i,2),'x','MarkerSize',6)
            hold on;
        end
        misLabel = labels~=yhat;
        plot(Y(misLabel,1),Y(misLabel,2),'ko','MarkerSize',8)
        if nClass==2
            legend('Class 1','Class 2','Class 1 Test','Class 2 Test','Misclassified')
        elseif nClass==3
            legend('Class 1','Class 2','Class 3',...
                        'Class 1 Test','Class 2 Test','Class 3 Test','Misclassified')
        end
        plotClassBoundary([X;Y],labels,@(X)naiveBayesPredict(X,model),false);
        set(gca,'FontSize',18)

        C = confusionmat(labels, yhat);
        plotConfMat(C);
        
    case 'nonparametric' %takes a long time for computation
        %generate random multivariate Gaussian data
        A1 = 0.2*randn(2,2);
        A2 = 0.15*rand(2,2);
        nClass1 = 400;
        nClass2 = 350;
        m1 = 3.4; m2 = 4;
        Y = mvnrnd(0.1*randn(1,2)+m1,(A1*A1'),nClass1);
        Y = [Y;mvnrnd(0.1*randn(1,2)+m2,(A2*A2'),nClass2)];
        labels = [ones(nClass1,1); 2*ones(nClass2,1)]; %create label vector
        if nClass > 2
            nClass3 = 375;
            m3 = 2.6;
            A3 = 0.3*rand(2,2);
            Y = [Y;mvnrnd(0.1*rand(1,2) + m3, A3*A3',nClass3)]; 
            labels = [ones(nClass1,1); 2*ones(nClass2,1);...
                        3*ones(nClass3,1)]; %create label vector
        end

        figure(1);
        for i=1:nClass
            plot(Y(labels==i,1),Y(labels==i,2),'.','MarkerSize',12)
            hold on;
        end

        model = naiveBayesFit(Y,labels,fcn);
        X = Y;

        %create different simulated data from same distributions
        %generate new multivariate Gaussian data for testing
        nClass1 = 200;
        nClass2 = 250;
        Y = mvnrnd(0.1*randn(1,2)+m1+0.05,(A1*A1'),nClass1);
        Y = [Y;mvnrnd(0.1*randn(1,2)+m2-0.05,(A2*A2'),nClass2)];
        labels = [ones(nClass1,1); 2*ones(nClass2,1)]; %create label vector
        if nClass > 2
            nClass3 = 375;
            Y = [Y;mvnrnd(0.1*rand(1,2) + m3-0.05, A3*A3',nClass3)]; 
            labels = [ones(nClass1,1); 2*ones(nClass2,1);...
                        3*ones(nClass3,1)]; %create label vector
        end
        tic
        [yhat,P] = naiveBayesPredict(Y,model);
        toc

        figure(1);
        for i=1:nClass
            plot(Y(labels==i,1),Y(labels==i,2),'x','MarkerSize',6)
            hold on;
        end
        misLabel = labels~=yhat;
        plot(Y(misLabel,1),Y(misLabel,2),'ko','MarkerSize',8)
        if nClass==2
            legend('Class 1','Class 2','Class 1 Test','Class 2 Test','Misclassified')
        elseif nClass==3
            legend('Class 1','Class 2','Class 3',...
                        'Class 1 Test','Class 2 Test','Class 3 Test','Misclassified')
        end
        plotClassBoundary([X;Y],labels,@(X)naiveBayesPredict(X,model),false);
        set(gca,'FontSize',18)

        C = confusionmat(labels, yhat);
        plotConfMat(C);
        
        
end


function plotConfMat(varargin)
%
%   usage: 
%   plotConfMat(confmat) plots the confmat with integers 1 to n as class labels
%   plotConfMat(confmat, labels) plots the confmat with the specified labels
%
%   Arguments
%   confmat:            a square confusion matrix
%   labels (optional):  vector of class labels

% number of arguments
switch (nargin)
    case 0
       confmat = 1;
       labels = {'1'};
    case 1
       confmat = varargin{1};
       labels = 1:size(confmat, 1);
    otherwise
       confmat = varargin{1};
       labels = varargin{2};
end

confmat(isnan(confmat))=0; % in case there are NaN elements
numlabels = size(confmat, 1); % number of labels

% calculate the percentage accuracies
confpercent = 100*confmat./repmat(sum(confmat, 1),numlabels,1);

% plotting the colors
figure,

imagesc(confpercent);
title(sprintf('Accuracy: %.2f%%', 100*trace(confmat)/sum(confmat(:))));
ylabel('Output Class'); xlabel('Target Class');

% set the colormap
colormap(flipud(gray));

% Create strings from the matrix values and remove spaces
textStrings = num2str([confpercent(:), confmat(:)], '%.1f%%\n%d\n');
textStrings = strtrim(cellstr(textStrings));

% Create x and y coordinates for the strings and plot them
[x,y] = meshgrid(1:numlabels);
hStrings = text(x(:),y(:),textStrings(:), ...
    'HorizontalAlignment','center');

% Get the middle value of the color range
midValue = mean(get(gca,'CLim'));

% Choose white or black for the text color of the strings so
% they can be easily seen over the background color
textColors = repmat(confpercent(:) > midValue,1,3);
set(hStrings,{'Color'},num2cell(textColors,2));

% Setting the axis labels
set(gca,'XTick',1:numlabels,...
    'XTickLabel',labels,...
    'YTick',1:numlabels,...
    'YTickLabel',labels,...
    'TickLength',[0 0]);
end





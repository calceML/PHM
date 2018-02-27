function model=RandomForestFit(X,y,NumTrees,NFeature)
% Fit a random forest for classification
%           model=RandomForestFit(X,y)
% Input:
%     X are independent variables
%     y is the response (label)
%     NumTrees is the number of trees. By default 20 trees are used.
%     NFeature is the number of features to select at random for each
%     decision split. Default value is sqrt(number of features)
% Output:
%     model is the fitted model.
%     class one data are plotted as blue dots; class two data are plotted
%     as red dots.

model=TreeBagger(NumTrees,X,y,'NumPredictorsToSample',NFeature);

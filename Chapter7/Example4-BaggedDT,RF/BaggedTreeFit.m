function model=BaggedTreeFit(X,y,NumTrees)
% Fit bagged trees for classification
%           model=BaggedTreeFit(X,y)
% Input:
%     X are independent variables
%     y is the response (label)
%     NumTrees is the number of trees. By default 20 trees are used.
% Output:
%     model is the fitted model.

model=TreeBagger(NumTrees,X,y,'NumPredictorsToSample','all');

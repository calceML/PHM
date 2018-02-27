function model=DecisionTreeFit(X,y)
% Fit a deicision tree for classification
%           model=DecisionTreeFit(X,y)
% Input:
%     X are independent variables
%     y is the response (label)
% Output:
%     model is the fitted model.

model=fitctree(X,y);




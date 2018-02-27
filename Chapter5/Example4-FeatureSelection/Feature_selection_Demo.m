% Feature_selection_Demo

% Load the data: 4000 features, 2 classes.
clear
load ovariancancer;
X=obs;
y=[ones(121,1);2*ones(216-121,1)];

% Perform feature selection
method={'ks','t','f','lda','svm'};
for k=1:length(method)
    ind=feature_selection(X,y,method{k});
    subplot(3,2,k)
    plot(ind,'.')
    title(method{k})
    xlabel('Feature Index');
    ylabel('Rank of Importance')
end

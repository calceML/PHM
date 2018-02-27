function ind=feature_selection(X,y,method)
% Select features for classification purpose.
%             ind=feature_selection(X,y,method)
% Input:
%     X are independent variables
%     y is the resposne (label)
% method:
%     'ks': Use two-sample KS test. It is non-parametricand does not assume distributions on the data.
%     't': Use two-sample t test. It assumes Gaussian distribution, and tests the separation by the mean.
%     'f': Use two-sample F test. It assumes Gaussian distribution, and tests the separation by the variance.
%     'lda': Use linear discriminant analysis (LDA), and assums Gaussian distribution.
%     'svm': Use support vector machine (SVM) with linear kernel, and does not assums distributions.
% Output:
%     ind is the rank of the importance of features feature.

switch method
    case 'ks'
        [pm,p_sort,ind]=feature_select_ks(X,y);
    case 't'
        [pm,p_sort,ind]=feature_select_t(X,y);
    case 'f'
        [pm,p_sort,ind]=feature_select_F(X,y);
    case 'lda'
        ind=feature_select_LDA(X,y);
    case 'svm'
        ind=feature_select_svm(X,y);
end

function [pm,p_sort,ind]=feature_select_ks(X,y)
% Select features for classification purpose using two-sample ks test.
% ks test is nonparametric so it does not have assumptions on the data.
%             pm=feature_select_ks(X,y)
% Input:
%     X are independent variables
%     y is the resposne (label)
% Output:
%     pm is the mean of p values of all pairs of classes.
%     Smaller pm value indicates better separation power of that feature.
%     p_sort is the sorted pm in ascending order.
%     ind is the corresponding feature index for p_sort.

label=unique(y);        % list label names
C=combnk(label,2);      % draw combinations of classes
T=size(C,1);            % Total number of combinations of classes
FE=size(X,2);           % Total number of features
P=zeros(T,FE);   
for k=1:T
    p=zeros(1,FE);
    for m=1:FE
        [h,p(1,m)]=kstest2(X(y==C(k,1),m),X(y==C(k,2),m));
    end
    P(k,:)=p;
end
pm=mean(P,1);
[p_sort,ind]=sort(pm);

function [pm,p_sort,ind]=feature_select_t(X,y)
% Select features for classification purpose using two-sample t test.
% Normality is assumed for each class. This function select features from
% their capability of differentiating different classes from mean values.
%             pm=feature_select_t(X,y)
% Input:
%     X are independent variables
%     y is the resposne (label)
% Output:
%     pm is the mean of p values of all pairs of classes.
%     Smaller pm value indicates better separation power of that feature.
%     p_sort is the sorted pm in ascending order.
%     ind is the corresponding feature index for p_sort.

label=unique(y);        % list label names
C=combnk(label,2);          % draw combinations of classes
T=size(C,1);            % Total number of combinations of classes
P=zeros(T,size(X,2));
for k=1:T
    [h,p,ci,stat]=ttest2(X(y==C(k,1),:),X(y==C(k,2),:));
    P(k,:)=p;
end
pm=mean(P,1);
[p_sort,ind]=sort(pm);

function [pm,p_sort,ind]=feature_select_F(X,y)
% Select features for classification purpose using two-sample F test.
% F test assumes the data follow Gaussian distribution. It tests if two
% samples have the sampe variance. Thus, this function focus on selecting
% features capable of separating the classes by differentiating the
% variances.
%             pm=feature_select_F(X,y)
% Input:
%     X are independent variables
%     y is the resposne (label)
% Output:
%     pm is the mean of p values of all pairs of classes.
%     Smaller pm value indicates better separation power of that feature.
%     p_sort is the sorted pm in ascending order.
%     ind is the corresponding feature index for p_sort.

label=unique(y);        % list label names
C=combnk(label,2);      % draw combinations of classes
T=size(C,1);            % Total number of combinations of classes
FE=size(X,2);           % Total number of features
P=zeros(T,FE);   
for k=1:T
    p=zeros(1,FE);
    for m=1:FE
        [h,p(1,m)]=vartest2(X(y==C(k,1),m),X(y==C(k,2),m));
    end
    P(k,:)=p;
end
pm=mean(P,1);
[p_sort,ind]=sort(pm);

function ind=feature_select_LDA(X,y)
% Select features for classification purpose using LDA.
%             ind=feature_select_svm(X,y)
% Input:
%     X are independent variables
%     y is the resposne (label)
% Output:
%     ind is the feature index ranked from top to bottom.

model=fitcdiscr(X,y);
[Delta_sort,ind]=sort(model.DeltaPredictor);

function ind=feature_select_svm(X,y)
% Select features for classification purpose using linear SVM.
%             ind=feature_select_svm(X,y)
% Input:
%     X are independent variables
%     y is the resposne (label)
% Output:
%     ind is the feature index ranked from top to bottom.

label=unique(y);        % list label names
C=combnk(label,2);          % draw combinations of classes
T=size(C,1);            % Total number of combinations of classes
for k=1:T
    yy=[y(y==C(k,1),:);y(y==C(k,2),:)];
    XX=[X(y==C(k,1),:);X(y==C(k,2),:)];
    Model=fitcsvm(XX,yy);
    [weight_sort,index]=sort(Model.Beta);
    w_index(:,k)=index;
end
w_index_sum=sum(w_index,2);
[weight_sum_sort,ind]=sort(w_index_sum);


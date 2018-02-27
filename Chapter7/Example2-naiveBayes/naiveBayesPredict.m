function [class,P] = naiveBayesPredict(X,model)
%Perform naive Bayes classification
%Inputs
%   X: n x d data matrix, d number of features, n number of observations
%   model: 
%       mu:     K x d matrix where each row is a class centroid
%       V:      K x d matrix where each column entry is the variance within a 
%                   dimension and rows are classes
%       g:      1 x K vector of priors
%       dist:   the distribution used


switch model.dist
    case 'gaussian'
        mu = model.mu;
        V = model.V;
        g = model.g;
        K = length(model.g); %number of classes
        [n,d] = size(X);
        P = zeros(n,K);
        for i=1:K
            Xc = bsxfun(@minus,X,mu(i,:)); %center data
            Xc = bsxfun(@times,Xc.^2,(1./V(i,:)).^2); %square to use in Norm dist.
            Xc = -0.5*(sum(Xc,2)); %each column entry is independent
            c = 1/sqrt(((2*pi)^(d/2)*prod(V(i,:)))); %normalization constant

            %column vector representing prob of each data point given class i
            P(:,i) = c*exp(Xc);
        end
        
    case 'lognormal'
        mu = model.mu;
        V = model.V;
        g = model.g;
        K = length(model.g); %number of classes
        [n,d] = size(X);
        P = zeros(n,K);
        x = log(X);
        for i=1:K
            c = 1/sqrt(((2*pi)^(d/2)*prod(V(i,:)))); 
            c = c./prod(X,2); %normalization constant
            Xc = bsxfun(@minus,x,mu(i,:)); %center data
            Xc = bsxfun(@times,Xc.^2,(1./V(i,:)).^2); %square to use in Norm dist.
            Xc = -0.5*(sum(Xc,2)); %each column entry is independent
            
            %column vector representing prob of each data point given class i
            P(:,i) = c.*exp(Xc);
        end
        
    case 'weibull'
        a = model.a;
        b = model.b;
        g = model.g;
        K = length(model.g); %number of classes
        n = size(X,1);
        P = zeros(n,K);
        for i=1:K
            %Perform analysis on log of weibull to avoid small numbers
                % log(b*a^(-b)) + (b-1)*log(x) - (x/a)^b
            c = log(b(i,:).*a(i,:).^(-b(i,:)));  %log coefficients
            c = sum(c); %sum log coefficients for normalization constant
            
            x1 = bsxfun(@times,(b(i,:)-1),log(X)); 
            x1 = sum(x1,2); %sum across all dimensions
            
            x2 = bsxfun(@times,X,1./a(i,:));
            x2 = bsxfun(@power,x2,b(i,:)); 
            x2 = sum(x2,2); %sum across all dimensions
            
            x = bsxfun(@plus,x1-x2,c); %log of the probabilities
            
            %column vector representing prob of each data point given class i
            P(:,i) = exp(x); %transform from logP to P
        end
        
    case 'nonparametric'
        g = model.g;
        [n,d] = size(X);
        Y = model.data(:,1:end-1);
        y = model.data(:,end);
        K = unique(y); %number of classes
        f = zeros(size(X));
        P = zeros(n,length(K));
        for i=K'
            for j=1:d
                f(:,j) = ksdensity(Y(y==K(i),j),X(:,j),...
                                    'kernel','epanechnikov');
            end
            P(:,i) = prod(f,2);
        end
        
        
end

P = bsxfun(@times,P,g); %multiply by prior
P = bsxfun(@times,P,1./sum(P,2)); %divide by sum to turn into proper probabilities

[~,class] = max(P,[],2); %take max probability

end

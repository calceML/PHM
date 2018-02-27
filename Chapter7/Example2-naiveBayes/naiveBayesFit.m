function model = naiveBayesFit(X,labels,distribution)

%build model
K = unique(labels); %classes
[n,D] = size(X); %total number of dimensions/features

switch distribution
    case 'gaussian'
        mu = zeros(length(K),D);
        V = zeros(length(K),D);
        g = zeros(1,length(K));
        for i=1:length(K)
            x = X(labels==K(i),:);
            nk = size(x,1); %number of obs in class i
            g(i) = nk/n; %prior for class i, can also set equal priors
            [mu(i,:),V(i,:)] = normfit(x); %get mean, variance for each dim. and ith class
        end
        model.mu = mu;
        model.V = V;
        model.g = g;
    
    case 'lognormal'
        mu = zeros(length(K),D);
        V = zeros(length(K),D);
        g = zeros(1,length(K));
        for i=1:length(K)
            x = X(labels==K(i),:);
            nk = size(x,1); %number of obs in class i
            g(i) = nk/n; %prior for class i, can also set equal priors
            for j=1:D
                a = lognfit(x(:,j)); %get mean, variance for each dim. and ith class
                mu(i,j) = a(1);
                V(i,j) = a(2);
            end
        end
        model.mu = mu;
        model.V = V;
        model.g = g;
        
    case 'weibull'
        a = zeros(length(K),D);
        b = zeros(length(K),D);
        g = zeros(1,length(K));
        for i=1:length(K)
            x = X(labels==K(i),:);
            nk = size(x,1); %number of obs in class i
            g(i) = nk/n; %prior for class i, can also set equal priors
            for j=1:D
                p = wblfit(x(:,j)); %get mean, variance for each dim. and ith class
                a(i,j) = p(1);
                b(i,j) = p(2);
            end
        end
        model.a = a;
        model.b = b;
        model.g = g;
        
    case 'nonparametric'
        g = zeros(1,length(K));
        for i=1:length(K)
            x = X(labels==K(i),:);
            nk = size(x,1); %number of obs in class i
            g(i) = nk/n; %prior probabilities
        end
        model.data = [X,labels];
        model.g = g;
end


model.dist = distribution;


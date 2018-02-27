function k = calcekfun(data, kernelmethod, parameter)

if ((nargin==2) & strcmp(kernelmethod, 'linear'))
    k = data*data';
elseif ((nargin==3) & strcmp(kernelmethod, 'polynomial'))
    k = data*data'+1;
    k = k.^parameter;
elseif ((nargin==3) & strcmp(kernelmethod, 'gaussian'))
    nobservations = size(data, 1);
    k = data*data'/parameter^2;
    d = diag(k);
    k = k - ones(nobservations,1)*d'/2;
    k = k - d*ones(1, nobservations)/2;
    k = exp(k);
end

end
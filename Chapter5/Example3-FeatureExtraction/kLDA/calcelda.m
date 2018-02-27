function ldc = calcelda(data, varargin)

% Linear Discriminant Analysis (LDA) for supervised dimensionality
% reduction can work binary and multiclass. It is also detects class label
% automatically and it returns linear discriminant components that
% can be used to multiply any new data for reduction.
% Input:
%   data -- M x N matrix, where M is the dimensioinality of
%           features/variables and N is the dimensionality of observations
%           classlabels.
% Output:
%   ldc -- M (dimension of features/variables) x R transformation matrix,
%          where R is the dimensionality of reduced space. The default
%          value of R is C-1, where C is the number of class labels.

p = inputParser;
p.CaseSensitive = false;    % Names are not sensitive to case: 'a' matches 'A'

defaultPlotXDim = 'LDC1';
defaultPlotYDim = 'LDC2';

addRequired(p, 'data', @ismatrix);
addParameter(p, 'plotxdim', defaultPlotXDim, @ischar);
addParameter(p, 'plotydim', defaultPlotYDim, @ischar);
parse(p, data, varargin{:});

inputdata = p.Results.data(:,1:end-1)';
datalabels = p.Results.data(:,end);
classtypes = unique(datalabels);
nclasses = length(classtypes);
[M, N] = size(inputdata);

%% LDA
classmeans = cell(1, nclasses);
classcovs = cell(1, nclasses);
ncobservations = zeros(1, nclasses);

% compute mean, covariance matrix, and number of observations for each class
for m=1 : nclasses
    classdata = inputdata(:, datalabels==classtypes(m));
    classmeans(m) = {mean(classdata, 2)};
    classcovs(m) = {cov(classdata')};
    ncobservations(m) = size(classdata, 2);
end

% compute mean of all the observatioins
globalmean = zeros(M, 1);
for m=1 : nclasses
    globalmean = globalmean+classmeans{m};
end
globalmean = globalmean/nclasses;

% compute between-class and within-class scatters, respectively
Sb = zeros(M, M);
Sw = zeros(M, M);

for m=1 : nclasses
    Sb = Sb + ncobservations(m).*(classmeans{m}-globalmean)*(classmeans{m}-globalmean)';
    Sw = Sw + classcovs{m};
end

% to reduce dimensionality, it is necessaryt to find "transmat" that
% maximize the ratio of between-class scatter and within-class scatter.
% "transmat" has to satisfy:
% (1) distance between the class means: the larger the better
% (2) the variance of each class: the smaller the better
% thus invSw is computed and the projection "transmat" that maximizes the
% ratio. This problem is converted to an Eigen vector problem for
% 1<=R<=C-1, where transmat = [projv1 | projv2 | ... |projvC-1]

invSw = inv(Sw);
[eigvec, eigval] = eig(invSw*Sb);
eigval = diag(eigval);

% sort eigval and select the top eigenvectors associated with the top
% eigenvalues as follows
[sorteigval, sorteigvalidx] = sort(eigval, 'descend');
ldc = eigvec(:, sorteigvalidx(1:nclasses-1));


plotxdim = str2num(p.Results.plotxdim(isstrprop(p.Results.plotxdim, 'digit')));
plotydim = str2num(p.Results.plotydim(isstrprop(p.Results.plotydim, 'digit')));
plotldcs = [plotxdim plotydim];

% plot projected data on linear discriminant components
calceldaplot(inputdata', datalabels, classtypes, ldc, plotldcs);

% save ldc
saveFolder = 'C:\CALCE-PHM\LDA';
if ~exist(saveFolder, 'dir');
    mkdir(saveFolder);
end

csvwrite([saveFolder '\ldc.csv'], pc);

end
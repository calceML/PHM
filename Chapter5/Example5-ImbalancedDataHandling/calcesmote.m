
clc, clear all, close all;

oversamplingmethod = 'SMOTE';

tmp = csvread('.\imbdata.csv');
data = tmp(:,1:end-1);
label = tmp(:,end);

dimx = 1;
dimy = 2;

% oversampling --
switch oversamplingmethod
    case 'SMOTE'
        uniqueclasses = unique(label);
        classidx = [1:length(uniqueclasses)];
        for m=1 : length(uniqueclasses)
            nsamplesperclass(m) = length(find(label==uniqueclasses(m)));
        end
        [~, majorityclassidx] = max(nsamplesperclass);
        
        minorityclassidx = [];
        cnt = 1;
        for m=1 : length(uniqueclasses)
            if classidx(m)~=majorityclassidx
                minorityclassidx(cnt) = classidx(m);
                cnt = cnt+1;
            end
        end
        
        % to show original data
        figure,
        grid on;
        color = ['c', 'r', 'b', 'k', 'm', 'y'];
        shape = ['o', 'h'];
        
        for m=1 : length(uniqueclasses)
            coloridx = mod(m,length(color))+1;
            shapeidx = 1;
            if (coloridx==1) && (shapeidx<=length(shape)) 
                shapeidx = shapeidx+1;
            end
               
            plotoption = [color(coloridx) shape(shapeidx)];
            
            hold on,
            plot(data(find(label==uniqueclasses(m)),dimx), data(find(label==uniqueclasses(m)),dimy), shape(shapeidx), 'MarkerFaceColor', color(coloridx)); 
            xlabel(['Dim: ' num2str(dimx)]), ylabel(['Dim: ' num2str(dimy)]);
            
            if m==majorityclassidx
                legendinfo{m} = ['Majority Class: ' num2str(uniqueclasses(m))];
            else
                legendinfo{m} = ['Minority Class: ' num2str(uniqueclasses(m))];
            end
        end
        legend(legendinfo);
        hold off;
        
        balancedata = [];
        balancelabel = [];
        majoritydata = data(find(label==uniqueclasses(majorityclassidx)),:);
        majoritylabel = label(find(label==uniqueclasses(majorityclassidx)));
        for m=1 : length(uniqueclasses)
            coloridx = mod(m,length(color))+1;
            shapeidx = 1;
            
            if m==majorityclassidx
                balancedata = [balancedata ; majoritydata];
                balancelabel = [balancelabel ; majoritylabel];
            else
                minoritydata = data(find(label==uniqueclasses(m)),:);
                minoritylabel = label(find(label==uniqueclasses(m)));
        
                nneighbors = 3;
                perc = floor((nsamplesperclass(majorityclassidx)/nsamplesperclass(m)-1)*100);
                smotedata = SMOTE(minoritydata, nneighbors, perc);
                smotelabel = uniqueclasses(m)*ones(length(smotedata),1);
                
                balancedata = [balancedata ; [minoritydata ; smotedata]];
                balancelabel = [balancelabel ; [minoritylabel ; smotelabel]];
                
                figure,
                grid on;
                plot(minoritydata(:,dimx), minoritydata(:,dimy), shape(shapeidx), 'MarkerFaceColor', color(coloridx));
                hold on, plot(smotedata(:,dimx), smotedata(:,dimy), 'sg', 'MarkerFaceColor', 'g');
                xlabel(['Dim: ' num2str(dimx)]), ylabel(['Dim: ' num2str(dimy)]);
                legend(['Minority Class: ' num2str(uniqueclasses(m))], [oversamplingmethod ': Synthetic Samples']);
                
            end
        end
        
    case 'ADASYN'
        uniqueclasses = unique(label);
        classidx = [1:length(uniqueclasses)];
        for m=1 : length(uniqueclasses)
            nsamplesperclass(m) = length(find(label==uniqueclasses(m)));
        end
        [~, majorityclassidx] = max(nsamplesperclass);
        
        minorityclassidx = [];
        cnt = 1;
        for m=1 : length(uniqueclasses)
            if classidx(m)~=majorityclassidx
                minorityclassidx(cnt) = classidx(m);
                cnt = cnt+1;
            end
        end
        
        % to show original data
        figure,
        grid on;
        color = ['c', 'r', 'b', 'k', 'm', 'y'];
        shape = ['o', 'h'];
        
        for m=1 : length(uniqueclasses)
            coloridx = mod(m,length(color))+1;
            shapeidx = 1;
            if (coloridx==1) && (shapeidx<=length(shape)) 
                shapeidx = shapeidx+1;
            end
               
            plotoption = [color(coloridx) shape(shapeidx)];
            
            hold on,
            plot(data(find(label==uniqueclasses(m)),dimx), data(find(label==uniqueclasses(m)),dimy), shape(shapeidx), 'MarkerFaceColor', color(coloridx)); 
            xlabel(['Dim: ' num2str(dimx)]), ylabel(['Dim: ' num2str(dimy)]);
            
            if m==majorityclassidx
                legendinfo{m} = ['Majority Class: ' num2str(uniqueclasses(m))];
            else
                legendinfo{m} = ['Minority Class: ' num2str(uniqueclasses(m))];
            end
        end
        legend(legendinfo);
        hold off;
        
        balancedata = [];
        balancelabel = [];
        majoritydata = data(find(label==uniqueclasses(majorityclassidx)),:);
        majoritylabel = label(find(label==uniqueclasses(majorityclassidx)));
        for m=1 : length(uniqueclasses)
            coloridx = mod(m,length(color))+1;
            shapeidx = 1;
            
            if m==majorityclassidx
                balancedata = [balancedata ; majoritydata];
                balancelabel = [balancelabel ; majoritylabel];
            else
                minoritydata = data(find(label==uniqueclasses(m)),:);
                minoritylabel = label(find(label==uniqueclasses(m)));
        
                nneighbors = 3;
                perc = floor((nsamplesperclass(majorityclassidx)/nsamplesperclass(m)-1)*100);
                smotedata = SMOTE(minoritydata, nneighbors, perc);
                smotelabel = uniqueclasses(m)*ones(length(smotedata),1);
                
                balancedata = [balancedata ; [minoritydata ; smotedata]];
                balancelabel = [balancelabel ; [minoritylabel ; smotelabel]];
                
                figure,
                grid on;
                plot(minoritydata(:,dimx), minoritydata(:,dimy), shape(shapeidx), 'MarkerFaceColor', color(coloridx));
                hold on, plot(smotedata(:,dimx), smotedata(:,dimy), 'sg', 'MarkerFaceColor', 'g');
                xlabel(['Dim: ' num2str(dimx)]), ylabel(['Dim: ' num2str(dimy)]);
                legend(['Minority Class: ' num2str(uniqueclasses(m))], [oversamplingmethod ': Synthetic Samples']);
                
            end
        
        
        
        
    case 'NONE'
end

%% SMOTE
function SMOTEd = SMOTE(data, k, N)
%   data = SMOTE(data, k, N)
%     @data - matrix of data. n_observations by n_features
%     @k - nearest neighbors
%     @model - smote precentage
%     
%   Run SMOTE on data.
% 
%   @Author: Gregory Ditzler (gregory.ditzler@gmail.com)

%     smote.m
%     Copyright (C) 2013 Gregory Ditzler
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


T = size(data, 1);
SMOTEd = [];

% If N is less than 100%, randomize the minority class samples as only a 
% random percent of them will be SMOTEd
if N < 100
    T = round((N/100)*T);
    N = 100;
end

for i = 1:T
    nnarray = [];
    synthetic = [];

    % determine the euclidean distance between the current minority sample
    % and reset of the other minority samples. then sort them in ascending
    % order
    for j = 1:T
        if i ~= j 
            euclid_dist(j,:) = data(i,:) - data(j,:);
        else
        % ignore the sample we are currently at from further calculations. 
        euclid_dist(j,:) = inf * ones(1, size(data,2));
        end
    end
    euclid_dist = sqrt(sum(euclid_dist.^2,2));
    euclid_dist2 = sort(euclid_dist,'ascend');

    % if we a really dealing with an imbalanced data set we may not have
    % enough samples to reduce to k nearest neighbors; instead grab all of
    % them 
    if length(euclid_dist2)<=k,
        knn = euclid_dist2;
        k = length(euclid_dist2);
    else
        knn = euclid_dist2(1:k);
    end

    % determine the k-nearest neighbors to the minority sample that we are
    % interested in.
    for j = 1:length(euclid_dist),
        if sum(euclid_dist(j)==knn)
            % the current distance in euclid_dist is a nearest neigbor of
            % the minority sample. so nnarray will have the indices of the
            % nearest neighbors in the minority instance array.
            nnarray(end+1) = j;
        end
    end

    % generate the synthetic samples
    newindex = 1; % keeps a count of number of synthetic samples generated
    N1 = round(N/100); % DO NOT OVERWRITE N!!!!

    while N1~=0
        % Choose a random number between 1 and k, call it nn. This step 
        % chooses one of the k nearest neighbors of i
        nn = round((k-1)*rand+1); % perform a linear conversion to scale the
                                  % nn paramter between 1 and k
        gap = rand;
        dif = data(nnarray(nn), :) - data(i, :);
        synthetic(newindex,:) = data(i, :) + gap*dif;

        newindex = newindex+1;
        N1 = N1-1;
    end
    
    SMOTEd = [SMOTEd; synthetic];
    clear euclid_dist euclid_dist2 N1 nnarray synthetic nnarray
end
end

%% ADASYN

function [ostrainset, ostrainlabel] = oversampling(trainset, trainlabel)

features0 = trainset(trainlabel==0,:);
features1 = trainset(trainlabel==1,:);
labels0 = trainlabel(trainlabel==0);
labels1 = trainlabel(trainlabel==1);

% ADASYN: set up ADASYN parameters and call the function:

adasyn_features                 = [features0; features1];
adasyn_labels                   = [labels0  ; labels1  ];
adasyn_beta                     = [];   %let ADASYN choose default
adasyn_kDensity                 = [];   %let ADASYN choose default
adasyn_kSMOTE                   = [];   %let ADASYN choose default
adasyn_featuresAreNormalized    = false;    %false lets ADASYN handle normalization
    
[adasyn_featuresSyn, adasyn_labelsSyn] = ADASYN(adasyn_features, adasyn_labels, adasyn_beta, adasyn_kDensity, adasyn_kSMOTE, adasyn_featuresAreNormalized);

ostrainset = [features0 ; features1 ; adasyn_featuresSyn];
ostrainlabel = [labels0 ; labels1 ; adasyn_labelsSyn];

% --



% PLOTTING:

%plot input data:
figure;
hold on;
plot(features0(:,1), features0(:,2), 'r.');
plot(features1(:,1), features1(:,2), 'b.');
title('input point sets');
xlabel('feature_1');
ylabel('feature_2');
axis('equal');
hold off;

%plot synthesized examples in green:
figure;
hold on;
plot(adasyn_featuresSyn(:,1), adasyn_featuresSyn(:,2), 'g.');
plot(features0(:,1), features0(:,2), 'r.');
plot(features1(:,1), features1(:,2), 'b.');
title('input point sets and synthetic points generated by ADASYN');
xlabel('feature_1');
ylabel('feature_2');
axis('equal');
hold off;

end


function [out_featuresSyn, out_labelsSyn] = ADASYN(in_features, in_labels, in_beta, in_kDensity, in_kSMOTE, in_featuresAreNormalized)
%this function implements the ADASYN method as proposed in the following
%paper:
%
%[1]: H. He, Y. Bai, E.A. Garcia, and S. Li, "ADASYN: Adaptive Synthetic
%Sampling Approach for Imbalanced Learning", Proc. Int'l. J. Conf. Neural
%Networks, pp. 1322--1328, (2008).
%
%the implementation follows the notation and equation numbers given in
%section 3.1.4 of another paper:
%
%[2]: H. He and E.A. Garcia, "Learning from imbalanced data",
%Knowledge and Data Engineering, IEEE Transactions on 21, no. 9,
%pp. 1263--1284, (2009).
%
%
%the purpose of the ADASYN method is to improve class balance towards
%equally-sized classes for a given input dataset. this is achieved by
%synthetically creating new examples from the minority class via linear
%interpolation between existing minority class samples. this approach is
%known as the SMOTE method, cf. section 3.1.3 in [2]. ADASYN is an
%extension of SMOTE, creating more examples in the vicinity of the boundary
%between the two classes, than in the interior of the minority class.
%cf. the supplied script demo_ADASYN for an example of this.
%
%
% INPUTS:
%----------
%in_features:
%(N \times P) matrix of numerical features. each row is one example, each
%column is one feature, hence there are N examples with P features each.
%
%in_labels:
%boolean N-vector of labels, defining the classes to which the examples in
%in_features belong.
%
%in_beta [default: 1]:
%desired level of balance, where 0 means that the size of the minority
%class will not be changed, and 1 means that the minority class will be
%ADASYNed to have (approximately, due to rounding) the same size as the
%majority class. any value of in_beta between 0 and 1 provides a compromise
%between these two extremes.
%note that in_beta IS NOT the resulting class ratio, but a percentage of
%how much class balance is improved in comparison to the given class
%balance! 0 means nothing is improved in comparison to the given class
%balance and 1 means class sizes are perfectly equalized (except for small
%rounding-related deviations).
%
%in_kDensity [default: 5]:
%k for kNN used in ADASYN density estimation, i.e. in calculation of the
%\Gamma_i values in eq. (4) of reference [2]. this is the kNN call that
%regards examples from both classes.
%
%in_kSMOTE [default: 5]: 
%k for kNN used in subsequent SMOTE-style synthesis of new examples.
%this is the kNN call that regards only examples from the minority class.
%cf. eq. (1) in reference [2].
%
%in_featuresAreNormalized [default: true]:
%boolean indicating whether the features (i.e. the different columns) in
%in_features are already normalized to the same scale or not.
%by default normalized features are assumed as the input, i.e. the user is
%expected to apply a normalization method of choice before passing the data
%to the ADASYN function.
%the practical difference in the two values of in_featuresAreNormalized is
%the following:
%true:  Euclidean distance is used in all kNN calls. for reasonable
%       results, in_features should already be normalized when calling
%       ADASYN().
%false: standardized Euclidean distance (type "doc knnsearch" into MATLAB
%       console and look for 'seuclidean' for an explanation) is used in
%       all kNN calls. any normalization already present in in_features is
%       ignored, and instead unit variance normalization is used.
%       however, this does NOT modify the data. the normalization is only
%       applied within knnsearch.
%
%
% OUTPUTS:
%----------
%out_featuresSyn, out_labelsSyn:
%features and labels of ONLY the synthetically created examples.
%note that each entry of out_labelsSyn is the label of the minority class
%since only examples of the minority class are created.
%concatenating [in_features out_featuresSyn] and [in_labels out_labelsSyn]
%gives a new example set with the desired class balance.
%
%
%
%
%-------------------------------------------------------------------------
% Version: 1.0
% Date: 2015-04-17
% Author: Dominic Siedhoff
%-------------------------------------------------------------------------
%
%
%-------------------------------------------------------------------------
% Copyright (c) 2015 Dominic Siedhoff
%-------------------------------------------------------------------------
%
% License: This software may be freely used, shared and modified. It must
%          not be sold. It is provided without any explicit or implicit
%          warranty of any kind. This license text must be included with
%          every copy made.
%
%-------------------------------------------------------------------------




if nargin < 3 || isempty(in_beta)
    in_beta = 1;
end

if nargin < 4 || isempty(in_kDensity)
    in_kDensity = 5;
end

if nargin < 5 || isempty(in_kSMOTE)
    in_kSMOTE = 5;
end

if nargin < 6 || isempty(in_featuresAreNormalized)
    in_featuresAreNormalized = true;
end


if in_beta == 0
    %nothing needs to be done because beta==0 is defined to mean that
    %current class ratio is kept:
    out_featuresSyn = [];
    out_labelsSyn   = [];
    return;
end


if ~all(in_labels==0 | in_labels==1)
    error('ADASYN: in_labels may contain only the values 0 and 1.');
end


numZeros = sum(in_labels==0);
numOnes  = sum(in_labels==1);

if numOnes == numZeros
    %nothing needs to be done because if classes are already balanced, then
    %for any in_beta, this is already the desired result.
    out_featuresSyn = [];
    out_labelsSyn   = [];
    return;
else
    if numZeros > numOnes
        majLabel = logical(0);
        minLabel = logical(1);
    else
        majLabel = logical(1);
        minLabel = logical(0);
    end
end

%rename:
S = in_features;
clear in_features;

%feature sets by class:
Smin = S(in_labels==minLabel,:);
Smaj = S(in_labels==majLabel,:);

%eq (3):
G = (size(Smaj,1) - size(Smin,1)) * in_beta;
G = round(G);

%handle boundary cases:
if size(Smin,1)==0
    warning('ADASYN: there were no examples of the minority class in the data. hence balancing is not possible. Returning empty matrices.');
    out_featuresSyn = [];
    out_labelsSyn   = [];
    return;
end

if size(Smin,1)==1
    warning('ADASYN: there was only one example of the minority class in the data. Hence returning G copies of that single example for balancing.');
    out_featuresSyn = repmat(Smin, [G 1]);
    out_labelsSyn   = logical(minLabel * ones([G 1]));
    return;
end


if in_featuresAreNormalized
    knnDistance = 'euclidean';
else
    %IMPORTANT: using 'seuclidean' as the distance measure means that
    %standardized Euclidean distance is used, i.e. the standard deviation
    %of the coordinates is automatically divided away. hence, using
    %'seuclidean' instead of 'euclidean' saves the effort of normalizing
    %the feature values by a Z-transformation (0mean,1var).
    %cf. documentation of input parameter in_featuresAreNormalized for more
    %information.
    knnDistance = 'seuclidean';
end

%kNN for density estimation:
idcs = knnsearch_nonflat(S,Smin, 'K',in_kDensity+1, 'Distance',knnDistance);
%note: why in_kDensity+1? because Smin is a subset of S and hence all
%points in Smin have a trivial nearest neighbor in S with distance 0.
%but that neighbor is not interesting because it's the point from Smin
%itself. hence remove it:
idcs = idcs(:,2:end);


%compute the \Gamma values (eq. (4) in reference [2]):
Gamma = zeros([size(Smin,1) 1]);
for cmi=1:size(Smin,1)  %cmi: current minority example index
    cNNs = idcs(cmi,:);             %current NearestNeighbors
    cNNsLabels = in_labels(cNNs);   %labels of cNNs:
    cNNsLabelsMaj = (cNNsLabels == majLabel);
    cDelta = nnz(cNNsLabelsMaj);    %the Delta_i of eq. (4) in reference [2]
    
    %write Gamma, not yet normalized:
    Gamma(cmi) = cDelta / in_kDensity;

end

%normalize Gamma to give a distribution function:
if sum(Gamma)==0
    %if there is no class overlap w.r.t. these knn settings, create a
    %uniform distribution:
    Gamma = 1/length(Gamma) * ones(size(Gamma));
else
    %create nonuniform distribution:
    Gamma = Gamma / sum(Gamma);     %this makes it exactly eq. (4) in reference [2]
end

%compute g_i (eq. (5) in reference [2]):
%these g_i are the numbers of synthetic examples to be generated from each
%example in Smin
g = round(Gamma * G);

if sum(g)==0
    warning('ADASYN: Classes are already well-balanced (i.e. sum(g)==0). Returning empty matrices.');
    out_featuresSyn = [];
    out_labelsSyn   = [];
    return;
end

%with this g known, call the ADASYN_SMOTE subroutine...:
out_featuresSyn = ADASYN_SMOTE(Smin,g,in_kSMOTE,knnDistance);
%...and generate the labels:
out_labelsSyn = logical(minLabel * ones([size(out_featuresSyn,1) 1]));

end


function Ssyn = ADASYN_SMOTE(Smin,g,k,knnDistance)
%subroutine implementing SMOTE algorithm as it is to be used by function
%ADASYN(). cf. section 3.1.3 in the following paper for details:
%
%[2]: H. He and E.A. Garcia, "Learning from imbalanced data",
%Knowledge and Data Engineering, IEEE Transactions on 21, no. 9,
%pp. 1263--1284, (2009).
%
%
% INPUTS:
%----------
%
%Smin:
%minority set from ADASYN()
%
%g:
%minority example synthesis counts as computed by ADASYN()
%
%k:
%number of neighbors to be regarded in kNN for SMOTE algorithm
%
%knnDistance:
%distance function used in kNN for SMOTE algorithm. depends on ADASYN's
%parameter in_featuresAreNormalized. please type "help ADASYN" into
%MATLAB's console for more information
%
%
% OUTPUTS:
%----------
%Ssyn: set of synthetic examples created from input set Smin by applying
%the SMOTE algorithm


%determine nearest neighbors:
idcs = knnsearch_nonflat(Smin,Smin, 'K',k+1, 'Distance',knnDistance);
%note: why k+1? because we search kNNs of Smin in Smin itself and hence all
%points in Smin have a trivial nearest neighbor in Smin with distance 0.
%but that neighbor is not interesting because it's the point from Smin
%itself. hence remove it:
idcs = idcs(:,2:end);

%initialize output and writing target as an empty matrix
Ssyn = zeros([0 size(Smin,2)]);

%for every minority example xi...
for cei=1:size(Smin,1)  %cei: current example index
    
    %current minority example:
    xi = Smin(cei,:);
    
    %number of synthetic examples to be created from xi:
    gi = g(cei);
    
    %allocate space for gi examples to be created from xi:
    xiSyn = zeros(gi, size(Smin,2));
    
    %...iterate over synthetic examples to be created from xi and
    %random partner from set of nearest neighbors:
    for csi=1:gi    %csi: current synthetic example index
        
        %get random partner example from nearest neighbors of xi:
        %neighbor index:
        nIdx = idcs(cei, randi(size(idcs,2)));
        %neighbor:
        xiHat = Smin(nIdx,:);
        
        %create synthetic example as according to eq. (1) in reference [2]:
        delta = rand(1);
        xSyn = xi + delta * (xiHat - xi);
        
        %write it to xiSyn:
        xiSyn(csi,:) = xSyn;
        
    end
    
    %append examples synthesized from xi to overall synthetic example set:
    Ssyn = [Ssyn; xiSyn];
    
end

end



function [IDX,D] = knnsearch_nonflat(X,Y, varargin)
%wraps knnsearch from MATLAB's statistics toolbox.
%knnsearch_nonflat executes knnsearch only on the dimensions with nonzero
%standard deviation, i.e. the flat dimensions are not passed on to
%knnsearch. this prevents pdist2 from producing the following warning in
%the context of standardized Euclidean distance ('seuclidean') in the
%presence of flat dimensions:
%"Warning: Some columns of S are zeros."
%if this warning occurs, pdist2 (and as a consequence knnsearch) gives only
%bad dummy results because the standardized Euclidean distance can not be
%computed properly in the presence of flat dimensions.
%using knnsearch_nonflat prevents this by filtering out flat dimensions.

nonflatX = std(X) ~= 0;
nonflatY = std(Y) ~= 0;

nonflat = nonflatX & nonflatY;

[IDX,D] = knnsearch(X(:,nonflat), Y(:,nonflat), varargin{:});


end







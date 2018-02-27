function [X_ind,idx,Cen,c]=cmeans_clustering(X,k)
% This function performs fuzzy c-means clustering and plot the result
% [X_ind,idx,Cen,c]=cmeans_clustering(X,k)
% Input:
%   X is the data to be clustered. Every column is a dimension.
%   k is the number of clustered.
% Output:
%   X_ind is the data with cluster indeces appended as the final columns
%   idx: cluster indeces as the membership probability
%   Cen: cluster centroids
%   c: clusters. Each c is a cell with all observations of a cluster. An
%   observation belongs to a cluster when its membership to that cluster is
%   larger than the rest.

[Cen,U]=fcm(X,k);
idx=U';
X_ind=[X,idx];
maxU=max(U);
for q=1:k
    c{q}=X(U(q,:)==maxU,:);
end

function [center, U, obj_fcn] = fcm(data, cluster_n, options)
%FCM Data set clustering using fuzzy c-means clustering.
%   Copyright 1994-2002 The MathWorks, Inc. 

if nargin ~= 2 & nargin ~= 3,
	error('Too many or too few input arguments!');
end

data_n = size(data, 1);
in_n = size(data, 2);

% Change the following to set default options
default_options = [2;	% exponent for the partition matrix U
		100;	% max. number of iteration
		1e-5;	% min. amount of improvement
		1];	% info display during iteration 

if nargin == 2,
	options = default_options;
else
	% If "options" is not fully specified, pad it with default values.
	if length(options) < 4,
		tmp = default_options;
		tmp(1:length(options)) = options;
		options = tmp;
	end
	% If some entries of "options" are nan's, replace them with defaults.
	nan_index = find(isnan(options)==1);
	options(nan_index) = default_options(nan_index);
	if options(1) <= 1,
		error('The exponent should be greater than 1!');
	end
end

expo = options(1);		% Exponent for U
max_iter = options(2);		% Max. iteration
min_impro = options(3);		% Min. improvement
display = options(4);		% Display info or not

obj_fcn = zeros(max_iter, 1);	% Array for objective function

U = initfcm(cluster_n, data_n);			% Initial fuzzy partition
% Main loop
for i = 1:max_iter,
	[U, center, obj_fcn(i)] = stepfcm(data, U, cluster_n, expo);
	if display, 
		fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
	end
	% check termination condition
	if i > 1,
		if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, break; end,
	end
end

iter_n = i;	% Actual number of iterations 
obj_fcn(iter_n+1:max_iter) = [];

function U = initfcm(cluster_n, data_n)
%INITFCM Generate initial fuzzy partition matrix for fuzzy c-means clustering.
%   Copyright 1994-2002 The MathWorks, Inc. 

U = rand(cluster_n, data_n);
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);

function [U_new, center, obj_fcn] = stepfcm(data, U, cluster_n, expo)
%STEPFCM One step in fuzzy c-mean clustering.
%   Copyright 1994-2014 The MathWorks, Inc. 

mf = U.^expo;       % MF matrix after exponential modification
center = mf*data./(sum(mf,2)*ones(1,size(data,2))); %new center
dist = distfcm(center, data);       % fill the distance matrix
obj_fcn = sum(sum((dist.^2).*mf));  % objective function
tmp = dist.^(-2/(expo-1));      % calculate new U, suppose expo != 1
U_new = tmp./(ones(cluster_n, 1)*sum(tmp));

function out = distfcm(center, data)
%DISTFCM Distance measure in fuzzy c-mean clustering.
%	Roger Jang, 11-22-94, 6-27-95.
%       Copyright 1994-2002 The MathWorks, Inc. 

out = zeros(size(center, 1), size(data, 1));

% fill the output matrix

if size(center, 2) > 1,
    for k = 1:size(center, 1),
	out(k, :) = sqrt(sum(((data-ones(size(data, 1), 1)*center(k, :)).^2)'));
    end
else	% 1-D data
    for k = 1:size(center, 1),
	out(k, :) = abs(center(k)-data)';
    end
end

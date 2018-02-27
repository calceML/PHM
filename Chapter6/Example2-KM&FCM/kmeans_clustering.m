function [X_ind,idx,Cen,c]=kmeans_clustering(X,k)
% This function performs original k-means clustering and plot the result
%   [X_ind,idx,Cen,c]=kmeans_clustering(X,k)
% Input:
%   X is the data to be clustered. Every column is a dimension.
%   k is the number of clustered.
% Output:
%   X_ind is the data with cluster indeces appended as the final column
%   idx: cluster indeces
%   Cen: cluster centroids
%   c: clusters. Each c is a cell with all observations of a cluster.

[idx,Cen]=kmeans(X,k);
X_ind=[X,idx];
for q=1:k
    c{q}=X(idx==q,:);
end

clear
% Make data
X = [randn(100,2)*0.75+ones(100,2);
    randn(100,2)*0.5-ones(100,2)];
% Perform k-means clustering
k=2;
[X_ind,idx,Cen,c]=cmeans_clustering(X,k);

% dimensions to plot
dim1=1;
dim2=2;
% plot data of each cluster
for q=1:k
    plot(c{q}(:,dim1),c{q}(:,dim2),'.');
    hold on
end
% Plot cluster centroids
plot(Cen(:,dim1),Cen(:,dim2),'kx','MarkerSize',15,'LineWidth',3)
xlabel('x')
ylabel('y')



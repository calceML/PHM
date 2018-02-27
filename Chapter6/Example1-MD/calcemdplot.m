function calcemdplot(mdrefdata, mddata, threshold)

concatenated = [mdrefdata ; mddata];

nrefdatapoints = length(mdrefdata);
ndatapoints = length(mddata);

figure, plot([1:nrefdatapoints], mdrefdata, 'bo',...
    [nrefdatapoints+1 : nrefdatapoints+ndatapoints], mddata, 'ro');
hold on, plot([1:nrefdatapoints+ndatapoints], threshold*ones(nrefdatapoints+ndatapoints, 1), 'k-');
xlabel('# of data points'), ylabel('Mahalanobis distance');
legend('Reference', 'Test', 'Threshold');

end
function calceocsvmtestplot(x, plotxdim, plotydim)

%% To load trained a ocsvm model
loadFolder = 'C:\CALCE-PHM\OCSVMAD';
load([loadFolder '\ocsvmmodel.mat']);

%%

if isempty(plotydim)
    
    [~, scorePred] = predict(ocsvmmodel, x(:, plotxdim));
    anomalyidx = find(scorePred < 0);

    axisline = [1 : size(x,1)];
    figure, plot(axisline, x, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    hold on,
    plot(anomalyidx, x(anomalyidx), 'k.', 'MarkerSize', 10, 'LineWidth', 2);
    title('{\bf Anomaly Detection via One-Class SVM}');
    xlabel('# of observatioins');
    ylabel('Dim1');
    legend('Observations','Anomalies');
        
elseif ~isempty(plotydim)
    
    [~, scorePred] = predict(ocsvmmodel, x(:, [plotxdim, plotydim]));
    anomalyidx = find(scorePred < 0);
    
    figure, plot(x(:,plotxdim), x(:,plotydim), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    hold on,
    plot(x(anomalyidx, plotxdim), x(anomalyidx, plotydim), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
    title('{\bf Anomaly Detection via One-Class SVM}');
    xlabel(['Dim' num2str(plotxdim)]);
    ylabel(['Dim' num2str(plotydim)]);
    legend('Observations','Anomalies')
    
end

end
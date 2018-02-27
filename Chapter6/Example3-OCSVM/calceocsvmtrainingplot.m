function calceocsvmtrainingplot(x, plotxdim, plotydim)

%% To load trained a ocsvm model
loadFolder = 'C:\CALCE-PHM\OCSVMAD';
load([loadFolder '\ocsvmmodel.mat']);

%%

svInd = ocsvmmodel.IsSupportVector;

if isempty(plotydim)
    
    axisline = [1 : size(x,1)];
    figure, plot(axisline, x, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    hold on,
    plot(svInd, x(svInd), 'k.', 'MarkerSize', 10, 'LineWidth', 2);
    title('{\bf Anomaly Detection via One-Class SVM}');
    xlabel('# of observatioins');
    ylabel('Dim1');
    legend('Observations','Support Vectors');
    
elseif ~isempty(plotydim)
    
    h = 0.02;   % mesh grid step size
    [x1, x2] = meshgrid(min(x(:,plotxdim)) : h : max(x(:,plotxdim)),...
    min(x(:,plotydim)) : h : max(x(:,plotydim)));
    [~, score] = predict(ocsvmmodel, [x1(:), x2(:)]);
    scoreGrid = reshape(score, size(x1,1), size(x2,2));
    
    figure, plot(x(:,plotxdim), x(:,plotydim),'k.', 'MarkerSize', 10, 'LineWidth', 2)
    hold on
    plot(x(svInd, plotxdim), x(svInd, plotydim), 'ro', 'MarkerSize', 10, 'LineWidth', 2)
    contour(x1, x2, scoreGrid)
    colorbar;
    title('{\bf Anomaly Detection via One-Class SVM}');
    xlabel(['Dim' num2str(plotxdim)]);
    ylabel(['Dim' num2str(plotydim)]);
    legend('Observations','Support Vectors')
    
end

end
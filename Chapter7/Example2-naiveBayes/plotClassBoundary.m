function [X1sparse,X2sparse] = plotClassBoundary(X,y,predictFcn,newFigure)

% Plot data and the classification boundaries induced by the specified 
% predictFcn. 
%
% X          - an n-by-2 data matrix
% Y          - an n-by-1 matrix of class labels
% predictFcn - a handle to a function of the form yhat = predictFcn(Xtest)
% Optional args:
% stipple    - [true] If true, use stippling, else use pastel translucent
%              shading
%
% Maximum 14 classes - although this is easy to change
%
% Examples
% plotDecisionBoundary(X, y, @(Xtest)logregPredict(model, Xtest));
% predictFcn = @(Xtest) logregPredict(model, kernelRbfSigma(Xtest, X, rbfScale)); 
% plotDecisionBoundary(X, y, predictFcn);
%%
    
    stipple = true;
    colors = plotColors();
    contourProps = {'LineWidth', 2, 'LineColor', 'k'};
    resolution = 350;
    
    
    
%     nclasses = nunique(y);
    nclasses = length(unique(y));
    range = dataWindow(X);
    [X1grid, X2grid, yhat] = gridPredict(range, resolution, predictFcn);
    if size(yhat,2)>1, yhat = yhat(:,1); end
    [X1sparse, X2sparse, yhatSparse] = gridPredict(range, resolution / 2.5, predictFcn);
    [nrows, ncols] = size(X1grid);
%     Y = canonizeLabels(Y);
    if newFigure, figure; hold on; end
    for c=1:nclasses
        if ~stipple && ~isOctave
            h = contourf(X1grid, X2grid, reshape(yhat, nrows, ncols), 1:nclasses);
            set(h,'LineWidth',2, 'FaceAlpha',0.1);
        else
            X1sparse = X1sparse(:); X2sparse = X2sparse(:);
            plot(X1sparse(yhatSparse==c), X2sparse(yhatSparse==c), '.', 'Color', colors{c}, 'MarkerSize', 0.05);
            yhatGrid = reshape(yhat, nrows, ncols);
            contour(X1grid, X2grid, yhatGrid, 1:nclasses, contourProps{:});
        end
     
    end
    axis(range);
    box on;
    axis tight
    
end

function [X1, X2, yhat] = gridPredict(range, resolution, predictFcn)
       X1range = linspace(range(1), range(2), resolution);
       X2range = linspace(range(3), range(4), resolution);
       [X1, X2] = meshgrid(X1range, X2range);
       [yhat] = predictFcn([X1(:), X2(:)]);
%        yhat = canonizeLabels(yhat);
end    

function [colors, colorMap] = plotColors()


    lightblue = [55 155 255] / 255;
    orange    = [255 128 0   ] / 255;
    green     = [0   255 64  ] / 255;
    magenta   = [255 0   128 ] / 255;
    olivegreen    = [132 199 71  ] / 255;
    %cyan      = [61  220 176 ] / 255;
    yellow2    = [215 215 0   ] / 255;
    red1    = [255 25 0   ] / 255;
    brown     = [128 64  0   ] / 255;
    blue      = [0   0   255 ] / 255;
    red      = [255   0   0 ] / 255;
    black      = [0   0   0 ] / 255;
    gray      = [128   128   128 ] / 255;
     
    colors = { lightblue, orange, green,  red1,  ...
      brown, blue,  black};
    %colors = repmat(colors, 1, 5);
    
    colorMap.lightblue = lightblue;
    colorMap.orange    = orange;
    colorMap.green     = green;
    %colorMap.cyan      = cyan;
    colorMap.yellow    = yellow2;
    colorMap.magenta   = magenta;
    colorMap.olivegreen    = olivegreen;
    colorMap.brown     = brown;
    colorMap.blue      = blue;
end


function window = dataWindow(X)
% Find appropriate axis coordinates for the n-by-2 data matrix X
% The output can be passed directly to the axis command.  

    
    assert(size(X, 2) == 2);
    minX1 = min(X(:, 1));
    maxX1 = max(X(:, 1));
    minX2 = min(X(:, 2));
    maxX2 = max(X(:, 2));
    dx1 = 0.15*(maxX1 - minX1);
    dx2 = 0.15*(maxX2 - minX2);
    window = [minX1 - dx1, maxX1 + dx1, minX2 - dx2, maxX2 + dx2];
end
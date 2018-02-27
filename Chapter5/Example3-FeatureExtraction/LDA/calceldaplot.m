function calceldaplot(data, datalabels, classtypes, ldc, plotldcs);
    nldcs = size(ldc, 2);
    if nldcs==1
        projdata = data*ldc;
        cmap = jet(size(classtypes,1));
        figure,
        grid on;
        for m=1 : length(classtypes)
            hold on, plot(projdata(classtypes(m)==datalabels, 1), zeros(length(projdata(classtypes(m)==datalabels, 1)),1),...
                'o', 'MarkerFaceColor', cmap(m,:), 'MarkerEdgeColor', 'k');
            xlabel('LDC1');
            ylim([-0.5 0.5]);
            h = legend();
            str = get(h, 'string');
            newlegstr = ['Class ' num2str(classtypes(m))];
            h =  legend([str newlegstr]);
        end
        hold off;
    elseif nldcs>=2
        projdata = data*ldc(:,plotldcs);
        cmap = jet(size(classtypes,1));
        figure,
        grid on;
         for m=1 : length(classtypes)
            hold on, plot(projdata(classtypes(m)==datalabels, 1), projdata(classtypes(m)==datalabels, 2),...
                'o', 'MarkerFaceColor', cmap(m,:), 'MarkerEdgeColor', 'k');
            xlabel(['LDC' num2str(plotldcs(1))]), yxlabel(['LDC' num2str(plotldcs(2))]);
            h = legend();
            str = get(h, 'string');
            newlegstr = ['Class ' num2str(classtypes(m))];
            h =  legend([str newlegstr]);
        end
        hold off;
    end
end

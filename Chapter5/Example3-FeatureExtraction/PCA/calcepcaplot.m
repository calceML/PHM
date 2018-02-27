function calcepcaplot(data, pc, plotpcs);
    npcs = size(pc, 2);
    if npcs==1
        projdata = data*pc;
        cmap = jet(size(projdata,1));
        figure,
        grid on;
        for m=1 : size(projdata,1)
            hold on, plot(projdata(m,1), zeros(length(projdata(m,1)),1), 'o', 'MarkerFaceColor', cmap(m,:), 'MarkerEdgeColor', 'k');
        end
        xlabel(['PC1']);
        ylim([-0.5 0.5]);
        hold off;
    elseif npcs>=2
        projdata = data*pc(:,plotpcs);
        cmap = jet(size(projdata,1));
        figure,
        grid on;
        for m=1 : size(projdata,1)
            hold on, plot(projdata(m,1), projdata(m,2), 'o', 'MarkerFaceColor', cmap(m,:), 'MarkerEdgeColor', 'k');
        end
        xlabel(['PC' num2str(plotpcs(1))]), ylabel(['PC' num2str(plotpcs(2))]);
        hold off;
    end
end

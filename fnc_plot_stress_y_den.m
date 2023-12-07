function fnc_plot_stress_y_den(x,y,axLim_factor,UnloadLoad_stress,extrm_idx_new,fig_name,xLab,yLab,Leg,out_dir,res)

% text_loc = 0.02 * axLim_factor*max([max(y)]); % the distance of annotations from y=0 is a fraction of ylim_max!
text_loc = mean([axLim_factor*min([y]) axLim_factor*max([y])]);
if exist("extrm_idx_new","var")  && size(extrm_idx_new, 2) == 1 && any(extrm_idx_new ~= 0)
    for i = 1 : numel(UnloadLoad_stress)
        x3(i,:) = [UnloadLoad_stress(i),UnloadLoad_stress(i)];
        y3(i,:) = [axLim_factor*min([y]),axLim_factor*max([y])];
    end
end
f = figure ('Name',fig_name,'Position',[100 100 500 375]);
set(f,'defaulttextinterpreter','latex');
plot(x,y,'--o','LineWidth',1,'MarkerSize',5,'Color',[0.85 0.47 0.32])
hold on
if exist("extrm_idx_new","var")  && size(extrm_idx_new, 2) == 1 && any(extrm_idx_new ~= 0)
    for i = 1 : numel(UnloadLoad_stress)
        plot(x3(i,:),y3(i,:),'--','LineWidth',0.8,'Color',[0.5 0.5 0.5])
        text(x3(i,1), text_loc, num2str(x3(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Clipping', 'on','FontSize',9,'Interpreter','latex')
    end
end
hold off
xlabel(xLab,'FontSize',10,'Color','k','Interpreter','latex')
ylabel(yLab,'FontSize',10,'Color','k','Interpreter','latex')
legend(Leg,'FontSize',9,'Location','northeast','Interpreter','latex')
xlim([axLim_factor*min(x) axLim_factor*max(x)])
ylim([axLim_factor*min([y]) axLim_factor*max([y])])
ax = gca;
set(ax,'TickLabelInterpreter','latex')
grid on

% Save the figure to the desired formats
% exportgraphics(gcf, fullfile(out_dir, [fig_name '.jpg']), 'Resolution', res);
exportgraphics(gcf, fullfile(out_dir, [fig_name '.png']), 'Resolution', res);
% exportgraphics(gcf, fullfile(out_dir, [fig_name '.tif']), 'Resolution', res);
% print(fullfile(out_dir, [fig_name '.svg']), '-dsvg', '-r300');
% print(fullfile(out_dir, [fig_name '.eps']), '-depsc', '-r300');

end
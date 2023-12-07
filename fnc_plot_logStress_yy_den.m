function fnc_plot_logStress_strain_den(x,y,x2,y2,axLim_factor,UnloadLoad_stress,fig_name,xLab,yLab,Leg1,Leg2,out_dir,res)

% text_loc = 0.02 * axLim_factor*max([max(y);max(y2)]); % the distance of annotations from y=0 is a fraction of ylim_max!
text_loc = mean([min([min(y);min(y2)]) axLim_factor*max([max(y);max(y2)])]);
for i = 1 : numel(UnloadLoad_stress)
    x3(i,:) = [UnloadLoad_stress(i),UnloadLoad_stress(i)];
    y3(i,:) = [axLim_factor*min([min(y);min(y2)]),axLim_factor*max([max(y);max(y2)])];
end
f = figure ('Name',fig_name,'Position',[100 100 500 375]);
set(f,'defaulttextinterpreter','latex');
% Plot the noisy data
semilogx(x,y,'--o','LineWidth',1,'MarkerSize',5,'Color',[0.85 0.47 0.32])
hold on
% Plot the denoised data
semilogx(x2,y2,'-','LineWidth',2,'MarkerSize',10,'Color',[0.1 0.25 0.89])
for i = 1 : numel(UnloadLoad_stress)
    semilogx(x3(i,:),y3(i,:),'--','LineWidth',0.8,'Color',[0.5 0.5 0.5])
    text(x3(i,1), text_loc, num2str(x3(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Clipping', 'on','FontSize',9,'Interpreter','latex')
end
hold off
xlabel(xLab,'FontSize',10,'Color','k','Interpreter','latex')
ylabel(yLab,'FontSize',10,'Color','k','Interpreter','latex')
legend(Leg1,Leg2,'FontSize',9,'Location','northeast','Interpreter','latex')
xlim([min(x) axLim_factor*max(x)])
ylim([min([min(y);min(y2)]) axLim_factor*max([max(y);max(y2)])])
set (gca, 'YDir','reverse')
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
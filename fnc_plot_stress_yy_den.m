function fnc_plot_stress_strain_den(x,y,x2,y2,axLim_factor,UnloadLoad_stress,fig_name,xLab,yLab,Leg1,Leg2)

f = figure ('Name',fig_name,'Position',[100 100 500 375]);
set(f,'defaulttextinterpreter','latex');
text_loc = 0.02 * axLim_factor*max([max(y);max(y2)]); % the distance of annotations from y=0 is a fraction of ylim_max!
for i = 1 : numel(UnloadLoad_stress)
    x3(i,:) = [UnloadLoad_stress(i),UnloadLoad_stress(i)];
    y3(i,:) = [axLim_factor*min([min(y);min(y2)]),axLim_factor*max([max(y);max(y2)])];
end
% Plot the noisy data
plot(x,y,'--o','LineWidth',1,'MarkerSize',5,'Color',[0.85 0.47 0.32])
hold on
% Plot the denoised data
plot(x2,y2,'-','LineWidth',2,'MarkerSize',10,'Color',[0.1 0.25 0.89])
for i = 1 : numel(UnloadLoad_stress)
    plot(x3(i,:),y3(i,:),'--','LineWidth',0.8,'Color',[1 0.84 0])
    text(x3(i,1), 0.3, num2str(x3(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Clipping', 'on','FontSize',9,'Interpreter','latex')
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
end
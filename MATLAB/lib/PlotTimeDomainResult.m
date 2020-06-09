function [rmse,max_rmse]=PlotTimeDomainResult(name,rk,yk,ek,t)

    figure('Name',name,'Position',[680 490 737 488]);
    subplot(2,1,1);
    plot(t,yk,':r','LineWidth',1.5);hold on;grid on;
    plot(t,rk,'--b','LineWidth',1.5);
    axis auto;
    ylabel('(mm)','FontSize',14);
    legend({'yk','rk'},'FontSize',14,'Location','northeastoutside');

    subplot(2,1,2);grid on;
    plot(t, ek(1:length(t)),'LineWidth',1.5);grid on;
    axis auto;
    legend({'ek'},'FontSize',14,'Location','northeastoutside');
    ylabel('(mm)','FontSize',14);
    xlabel('time(s)','FontSize',16);

    rmse = rms(ek)*1000;
    max_rmse = max(ek)*1000;
    
    title_msg = sprintf("rmse: %.4f um", rmse);
    title(title_msg,'FontSize',16);

end
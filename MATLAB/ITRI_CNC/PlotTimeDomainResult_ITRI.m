function [rmse,max_rmse]=PlotTimeDomainResult_ITRI(name,rk,yk,ek,t,plant)

    figure('Name',name,'Position',[680 490 737 488]);
    subplot(2,1,1);
    plot(t,yk,':r','LineWidth',1.5);hold on;grid on;
    plot(t,rk,'--b','LineWidth',1.5);
    axis auto;
    ylabel('(count)','FontSize',14);
    legend({'yk','rk'},'FontSize',14,'Location','northeastoutside');
    title(name,'FontSize',16);
    
    subplot(2,1,2);grid on;
    plot(t, ek(1:length(t)),'LineWidth',1.5);grid on;
    axis auto;
    legend({'ek'},'FontSize',14,'Location','northeastoutside');
    ylabel('(count)','FontSize',14);
    xlabel('time(s)','FontSize',16);

    rmse = rms(ek);
    max_rmse = max(ek);
    
    title_msg = sprintf("rmse: %f counts (%.4f um)", rmse, rmse*plant.count2round*plant.round2mm*1000);
    title(title_msg,'FontSize',16);

end
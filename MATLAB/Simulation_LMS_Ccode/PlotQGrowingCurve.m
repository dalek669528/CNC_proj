function PlotQGrowingCurve(name,qk,Lq,t)

    colormsg = {[1 0 0];
                [0 1 0];
                [0 0 1];
                [1 1 0];
                [0 1 1];
                [1 0 1];
                [1 0.5 0.5];
                [0.5 1 0.5];
                [0.5 0.5 1];
                [1 1 0.5];
                [0.5 1 1];
                [1 0.5 1];
                [1 0 0.5];
                [1 0.5 0];
                [0.5 1 0];
                [0 1 0.5]};

    msg = cell(1,Lq);
    for i = 1:Lq
        msg(i) = {"q"+i};
    end

    figure('Name',name,'Position',[680 490 837 488]);

    subplot(1,3,[1,2]);
    h0 = plot(t,qk,'LineWidth',1);grid on;
    set(h0,{'Color'},colormsg(1:Lq));
    axis auto;
    xlabel('time(s)','FontSize',14);    
    title(name,'FontSize',16);

    legend(msg,'FontSize',14,'Location','northwestoutside');
    
    segment = 1:1000;
    subplot(1,3,3);
    h1 = plot(t(segment),qk(segment,:),'LineWidth',1);grid on;
    set(h1,{'Color'},colormsg(1:Lq));
    axis auto;
    xlabel('time(s)','FontSize',14);    
    title("1 sec",'FontSize',16);

    
end
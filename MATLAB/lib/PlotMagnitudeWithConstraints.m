function [mag]=PlotMagnitudeWithConstraints(name,sys,F)
    fmin = 0.01*2*pi;
    fmax = 500*2*pi;

    [mag,~,wout] = bode(sys,fmin:0.01:fmax);
    [Gm,Pm,Wcg,Wcp] = margin(sys);
    [Nc,~] = size(F);
    
    if Nc ~= 0
        Fpoint = F(:,1);
        Mpoint = 20*log10(F(:,2));
    end

    figure('Name',name+"-Bode Diagram");
    semilogx(wout,20*log10(mag(:,:)),'b');
    grid on;hold on;
    
    if Nc ~= 0
        area(Fpoint,Mpoint,'FaceColor',[1 0.9 1],'EdgeColor',[1 0.9 1]);
        plot(Fpoint,Mpoint,'-.*','Color',[1 0.4 0.4],'LineWidth',1);        
    end
    
    plot([fmin,fmax],[0,0],'black','LineWidth',0.5);
    
    semilogx(wout,20*log10(mag(:,:)),'b','LineWidth',1.5);
    
    axis([fmin fmax -inf inf]);
    axis('auto y');
    s1 = sprintf("GM= %.4f dB   (at %.4f rad/s)", 20*log10(Gm), Wcg);
    s2 = sprintf("Pm= %.4f deg (at %.4f rad/s)", Pm, Wcp);
    title({s1,s2});
    ylabel("Magnitude(dB)");
    xlabel("Frequency(rad/s)");
    
end
function fit = CheckFitnessFun(q, F)
    count = zeros(F.Nc,1);
    for i = 1:F.Nc
        alpha = F.alpha(i);
        beta = F.beta(:,i);
        Omega = F.Omega(:,:,i);
        gamma = F.gamma(i);
        g = alpha* q.'*Omega*q + beta.'*q + gamma;
%         if g>0
            count(i) = g;
%         end
    end
    fit = count;
end
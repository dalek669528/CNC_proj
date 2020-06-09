function Saveconstraints(type,F)
Lq = length(F.Omega(:,:,1));
fid_p=fopen(strcat(type,'_DesignedF.txt'),'w');
for i = 1:F.Nc
    beta = F.beta(:,i);
    Omega = F.Omega(:,:,i);
    Ochain = [Omega(Lq,1:length(Omega)-1).'; Omega(:,1)];
    gamma = F.gamma(i);
    weight = F.weight(i);

    fprintf(fid_p,'%.19f\t', beta);

    fprintf(fid_p,'%.19f\t', gamma);

    fprintf(fid_p,'%.19f\t', Ochain);
    
    fprintf(fid_p,'%.19f\t', weight);  
    if i ~= F.Nc
        fprintf(fid_p,'\r\n');
    end
end
fclose(fid_p);
end

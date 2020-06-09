%% Linear Fractional Transformation (LFT) function
function plantnew = LinearFractionalTransformation(plant)

    % NOTICE !!!
    % There is one minus before the model
    G = minreal(-1*canon(plant.v2p));
    pole = plant.assignedpole;
    F = place(G.A, -1*G.B, pole);
    H = place(G.A', -1*G.C', pole);
    L = H.'; 

    JA = G.A + G.B*F + L*G.C +L*G.D*F;
    JB = [-L G.B+L*G.D] ;
    JC = 1*[F; -(G.C+G.D*F)];
    JD = 1*[0 eye(length(G.D)); eye(length(G.D)) -1*G.D];    
    Jy = minreal(ss(JA,JB,JC,JD,plant.Ts));
    plant.Jy = Jy;

    plantnew = plant;
end

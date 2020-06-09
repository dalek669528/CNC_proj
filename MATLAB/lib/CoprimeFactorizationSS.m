%%  Youla Parameterization via Coprime Facterization 
%   input: id model & pole assignment
%   output: modified model
function newplant = CoprimeFactorizationSS(plant, gain)
    Ts = plant.Ts;
    G = minreal(canon(plant.v2p));
        % canon() : transfer the v2p to canonical state-space form
        % minreal() : do pole-zero cancellation
    pole = zeros(length(eig(G.A)),1);
        % length() : if a matrix A is m by n, length(A) will return max(m, n)
        % zeros(m, n) : create an matrix m by n, and all of the elements are zeros
    pole = complex(pole,0);
    eigA = eig(G.A);

    % z domain , all pole should be inside the unit circle
    for i = 1:length(G.A)
        if abs(eigA(i))>0.99
            pole(i) = eigA(i)*gain;
        else
            pole(i) = eigA(i);
        end
    end

    % reference to the coprime factorization
    F = place(G.A, -1*G.B, pole);
    H = place(G.A', -1*G.C', pole);
   
    M = minreal(ss(G.A+G.B*F, G.B, F, 1, Ts));
    N = minreal(ss(G.A+G.B*F, G.B, G.C+G.D*F, G.D, Ts));
    X = minreal(ss(G.A+H.'*G.C, H.', F, 0, Ts));
    W = minreal(ss(G.A+H.'*G.C, -G.B-H.'*G.D, F, 1, Ts));
    Bezout_identity = minreal(N*X+M*W); %should be 1
    
    plant.assignedpole = pole;
    plant.M = M;
    plant.N = N;
    plant.X = X;
    plant.W = W;
    plant.Bezout_identity = Bezout_identity;

    newplant = plant;
end

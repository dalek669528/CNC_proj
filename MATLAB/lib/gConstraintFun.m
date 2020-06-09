function allFrequencyconstraints = gConstraintFun(DesignedF,L,plant)
    [Nc,~] = size(DesignedF);
    alpha = zeros(1,Nc);
    beta = zeros(L,Nc);
    Omega = zeros(L,L,Nc);
    gamma = zeros(1,Nc);
    weight = zeros(1,Nc);
    N = plant.N;
    M = plant.M;
    X = plant.X;
    W = plant.W;
    P = plant.v2p;
    Ts = P.Ts;
    
    for i = 1:Nc
        % assign the weight for each constraint
        weight(i) = 0.0000005;
        w0 = 100000;

        % the frequency(rad/sec) of constraint
        wi = DesignedF(i,1);

        % the gain(unit) of constraint
        %     ki > 1    ->  lower bounded case
        % 0 < ki < 1    ->  upper bounded case
        ki = DesignedF(i,2);
        
        % get the value of each model at specific frequency
        % assign frequency(rad/s) and back gain(unit) and phase(deg)
        [magP,~] = bode(P,wi);
        [magM,phaseM] = bode(M,wi);
        [magN,phaseN] = bode(N,wi);
        [magX,phaseX] = bode(X,wi);
        [magW,phaseW] = bode(W,wi);
        valueM = magM*cosd(phaseM) + 1i*magM*sind(phaseM);
        valueN = magN*cosd(phaseN) + 1i*magN*sind(phaseN);
        valueX = magX*cosd(phaseX) + 1i*magX*sind(phaseX);
        valueW = magW*cosd(phaseW) + 1i*magW*sind(phaseW);
    
        % compute the all parameters
        % (alpha) * q.' * (Omega) * q + (beta).' * q + (gamma);
        rho = zeros(L,1);
        rho = complex(rho,0); % rho = complex(a,b) create a complex output, rho, from two real inputs, such that rho = a + bi
        omevector = zeros(L,1);
        % discrete frequency = continuos frequency * sampling time
        % range [0 , pi]
        for n = 0:L-1
            rho(n+1) = exp(-n*wi*1i*Ts);
            omevector(n+1) = cos(n*wi*Ts);
        end
        Omega(:,:,i) = toeplitz(omevector);
%         example of using toeplitz(r), assume r is a real vector :
%           r = [1 2 3];
%           toeplitz(r) = 1 2 3
%                         2 1 2
%                         3 2 1

        b = real(conj(valueW) * valueN * rho);
        c = real(conj(valueX) * valueM * rho);
        
        alpha(:,i) = 1; % useless
        beta(:,i) = (-2*(ki.^2*b+magP.^2*c))/((ki.^2-1)*(magN.^2)); 
        gamma(:,i) = ((ki.^2)*(magW.^2)-(magP.^2)*(magX.^2))/((ki.^2-1)*(magN.^2));

    end

    allFrequencyconstraints = struct('alpha',alpha,...
                                     'beta',beta,...
                                     'Omega',Omega,...
                                     'Nc',Nc,...
                                     'gamma',gamma,...
                                     'weight',weight,...
                                     'w0',w0);
end
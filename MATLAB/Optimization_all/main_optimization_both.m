clc;clear;close all;

%% Load the ID model from lib
ID_Model;
Designedplant = plantZ;

%% Design reference input
Reference_time = 15; % (sec)
Reference_amplitude = 50; % (mm)
Reference_frequency = 2*pi*1/15;
chirp_min = 0.01; % the lowest frequency
chirp_max = 1.06; % the highest frequency
inputdata = myChirpsine(Reference_time,Ts,chirp_min,chirp_max);
Reference_inputdata = inputdata;

%% Design the Frequency constraints
% new design (0.4952/5) Z       
DesignedF = [
            0.1     1000;
            1       100;
            10      10;
            100     1.1;
            300     0.8;
            400     0.56;
            500     0.39;
            600     0.26;
            700     0.16;
            800     0.09;
            900     0.07;
            1000    0.08;
            1300    0.1;        
            1500    0.09;
            1800    0.06; 
            2000    0.03;
            2500    0.07;
            3000    0.1            
            ];

[Ncc,~] = size(DesignedF);

%% Determine the length and pole
% all_pole = (25:5:75)/100;
% all_pole = [0.25 0.75];
% all_length = 5:1:20;
% all_length = [5 15 50 80];
all_pole = [0.4952];
all_length = [15];

%% Solve the QCQP problem
allfitness = zeros(Ncc,length(all_length),length(all_pole));
allQ = zeros(max(all_length),length(all_length),length(all_pole));
allstatus = cell(length(all_length),length(all_pole)); % cell array which can store any type of data
allrms = zeros(2,length(all_length),length(all_pole));
allmax = zeros(2,length(all_length),length(all_pole));
count_pole = 0;

for pole = all_pole
    count_pole = count_pole + 1;
% Do Coprime Factorization from lib
    plant = CoprimeFactorizationSS(Designedplant,pole);
    plant = LinearFractionalTransformation(plant);
    sim('Optimization_getVZ'); % do Simulink.
    
    count_L = 0;
    for L = all_length % L -> length of q
        count_L = count_L + 1;
% Time-domain coefficient %
        Tphi = zeros(L,1);
        Tcoe2_t = zeros(L,L);
        Tcoe1_t = zeros(1,L);
        Tcoe0_t = 0;
        for k = 1:length(vk)
            Tphi = [vk(k) ;Tphi(1:L-1)];
            Tcoe2_t = Tcoe2_t + Tphi*Tphi.';
            Tcoe1_t = Tcoe1_t + zk(k)*Tphi.';
            Tcoe0_t = Tcoe0_t + zk(k)^2; % sigma( zk(k)^2 )
        end
    
% Frequency domain coefficient %
        allFC = gConstraintFun(DesignedF,L,plant);
%  use cvx to solve qcqp %
        tic;
        cvx_begin 
        cvx_precision 
            variable q(L)
            minimize((1.0/2.0) * quad_form(q,(2.0)*Tcoe2_t) + (-2*Tcoe1_t)*q + Tcoe0_t)
            subject to
            for fi = 1:allFC.Nc
                (1.0/2.0) * quad_form(q,(2.0)*allFC.Omega(:,:,fi)) + allFC.beta(:,fi).'*q + allFC.gamma(fi) <= 0;
               % (1.0/2.0) * 10^(-5) * quad_form(q,(2.0)*allFC.Omega(:,:,fi)) + 10^(-5) * allFC.beta(:,fi).'*q + 10^(-5) * allFC.gamma(fi) <= 0;
            end
        cvx_end
        ggg = toc
        
        fprintf("pole = %.3f, L = %d, status: %s\n", pole, L,cvx_status);
        
% save data to txt %
        numberOfConstraints = allFC.Nc;
%         save("C:\Users\LAB_JHU\Desktop\ken\online solver for QCQP test\Totsu-rust_v0.4.0\solver\src\data\numberOfConstraints_.txt", 'numberOfConstraints', '-ascii');
        Q_length = length(q);
%         save("C:\Users\LAB_JHU\Desktop\ken\online solver for QCQP test\Totsu-rust_v0.4.0\solver\src\data\Q_length_.txt", 'Q_length', '-ascii');
        P_matrix = zeros(L, L);
        q_vector = zeros(L, 1);
        r_scalar = zeros(1, 1);
        for fi = 0:1:allFC.Nc
            if(fi == 0)
                P_matrix = Tcoe2_t * (2.0);
                q_vector = Tcoe1_t' * (-2.0);
                r_scalar = Tcoe0_t;
            else
                P_matrix = 10^(-5) * allFC.Omega(:,:,fi) * (2.0);
                q_vector = 10^(-5) * allFC.beta(:,fi);
                r_scalar = 10^(-5) * allFC.gamma(fi);
            end
%             save("C:\Users\LAB_JHU\Desktop\ken\online solver for QCQP test\Totsu-rust_v0.4.0\solver\src\data\P_matrix\P_matrix_" + fi + "_.txt", 'P_matrix', '-ascii');
%             save("C:\Users\LAB_JHU\Desktop\ken\online solver for QCQP test\Totsu-rust_v0.4.0\solver\src\data\q_vector\q_vector_" + fi + "_.txt", 'q_vector', '-ascii');
%             save("C:\Users\LAB_JHU\Desktop\ken\online solver for QCQP test\Totsu-rust_v0.4.0\solver\src\data\r_scalar\r_scalar_" + fi + "_.txt", 'r_scalar', '-ascii');
        end
        
% record data %
        allfitness(:,count_L,count_pole) = CheckFitnessFun(q,allFC);
        allstatus(count_L,count_pole) = {cvx_status};
        allQ(1:L,count_L,count_pole) = q;

% perfomance test %
        % q = [4670.21 -2916.04 2493.59 -737.039 1800.59 1362.52 511.299 967.203 1602.43 366.933 -471.903 1470.56 -2359.58 2397.5 -3554.25].';
        [ccA,ccB,ccC,ccD] = LFTExpandedSS(plant.Jy,q);
        Lq = L; % for simulink
        Qnum = q;
        
        for f=1
            Test_time = 30;
            Test_amplitude = 50;
            Test_frequency = 2*pi*f/15;
%             inputdata = myCosine(Test_time,Ts,Test_frequency,Test_amplitude);
            inputdata = myDesignedCosine(Test_time,Ts,Test_frequency,Test_amplitude);
            sim('Block_fixedQ_Ccode_ideal');
            allrms(f,count_L,count_pole) = rms(ek);
            allmax(f,count_L,count_pole) = max(abs(ek));
        end
        Test_inputdata = inputdata;
    end
end
 clear allFC f fi k L Lq pole q Tcoe0_t Tcoe1_t Tcoe2_t Tphi tout

%%
% save("Z_chirp1hz_lengthtest.mat",'Designedplant','all_length','all_pole','allfitness','allmax','allQ','allrms','allstatus','DesignedF','Reference_inputdata','Test_inputdata');

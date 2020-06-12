clc;clear;close all;

date = "0612";
if ~exist("C:/Users/dalek669528/Desktop/CNC_proj/MATLAB/data/" + date, 'dir')
   mkdir("C:/Users/dalek669528/Desktop/CNC_proj/MATLAB/data/" + date)
end

Q_length = 15;
poles = 0.301:0.001:0.6;
% DesignT_change_times = 1;
randn_rate = 0.02;
% cvx_lims = [9e-2, 6e-2, 3e-2, 1e-2, 8e-3, 5e-3, 3e-3, 5e-4, 1e-4, 5e-5, 1e-5, 1e-6];

iter_times = length(poles)*length(Q_length);
checkpoint = iter_times;
max_err_X = zeros(1, iter_times);
max_err_Z = zeros(1, iter_times);
max_mid_err_X = zeros(1, iter_times);
max_mid_err_Z = zeros(1, iter_times);
rmse_err_X = zeros(1, iter_times);
rmse_err_Z = zeros(1, iter_times);
rmse_mid_err_X = zeros(1, iter_times);
rmse_mid_err_Z = zeros(1, iter_times);

 
ID_Model_4;

%% Design reference input
Reference_time = 15; % (sec)
Reference_amplitude = 50; % (mm)
Reference_frequency = 2*pi*1/15;
chirp_min = 0.01; % the lowest frequency
chirp_max = 1.06; % the highest frequency
inputdata = myChirpsine(Reference_time,Ts,chirp_min,chirp_max);
Reference_inputdata = inputdata;


%% Design the Frequency constraints
DesignedF_ori = [
    0.1     1000;
    1       100;
    10      10;
    100     1.3;
    350     0.831;
    500     0.5;
    700     0.15;
    1000    0.1;
    1400    0.1;
    1800    0.1;
    2150    0.09;
    2500    0.08;
    3000    0.083
    ];
DesignedF = DesignedF_ori;
[Ncc,~] = size(DesignedF);

idx = 1;
% for cvx_lim = cvx_lims
for Designedplant = [plantZ, plantX, plantZ_2, plantX_2]
    %% Solve the QCQP problem
    allfitness = zeros(Ncc,length(Q_length),length(poles));
    allQ = zeros(max(Q_length),length(Q_length),length(poles));
    allstatus = cell(length(Q_length),length(poles)); % cell array which can store any type of data
    allrms = zeros(length(Q_length),length(poles));
    allmax = zeros(length(Q_length),length(poles));
    count_pole = 0;
    for pole = poles
        count_pole = count_pole + 1;
        % Do Coprime Factorization from lib
        plant = CoprimeFactorizationSS(Designedplant,pole);
        plant = LinearFractionalTransformation(plant);
        sim('Optimization_getVZ'); % do Simulink.

        count_L = 0;
        for L = Q_length % L -> length of q
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
            % use cvx to solve qcqp %
            tic;
            cvx_begin quiet
            % args = [ cvx_lim, cvx_lim*5, cvx_lim*10];
            % cvx_precision (args)
            cvx_precision low
                variable q(L)
                minimize((1.0/2.0) * quad_form(q,(2.0)*Tcoe2_t) + (-2*Tcoe1_t)*q + Tcoe0_t)
                subject to
                for fi = 1:allFC.Nc
                   (1.0/2.0) * quad_form(q,(2.0)*allFC.Omega(:,:,fi)) + allFC.beta(:,fi).'*q + allFC.gamma(fi) <= 0;
                    % (1.0/2.0) * 10^(-5) * quad_form(q,(2.0)*allFC.Omega(:,:,fi)) + 10^(-5) * allFC.beta(:,fi).'*q + 10^(-5) * allFC.gamma(fi) <= 0;
                end
            cvx_end
            ggg = toc;

            fprintf("pole = %.3f, L = %d, status: %s\n", pole, L, cvx_status);
            % save data to txt %
            numberOfConstraints = allFC.Nc;
            % save(".\data\numberOfConstraints_.txt", 'numberOfConstraints', '-ascii');
            % save(".\data\Q_length_.txt", 'Q_length', '-ascii');
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
                % save(".\data\P_matrix\P_matrix_" + fi + "_.txt", 'P_matrix', '-ascii');
                % save(".\data\q_vector\q_vector_" + fi + "_.txt", 'q_vector', '-ascii');
                % save(".\data\r_scalar\r_scalar_" + fi + "_.txt", 'r_scalar', '-ascii');
            end

            % record data %
            allfitness(:,count_L,count_pole) = CheckFitnessFun(q,allFC);
            allstatus(count_L,count_pole) = {cvx_status};
            allQ(1:L,count_L,count_pole) = q;

            % perfomance test %
            [ccA,ccB,ccC,ccD] = LFTExpandedSS(plant.Jy,q);
            Lq = L; % for simulink
            Qnum = q;
            Test_time = 30;
            Test_amplitude = 50;
            Test_frequency = 2*pi/15;
            inputdata = myCosine(Test_time,Ts,Test_frequency,Test_amplitude);
            % if isequal(Designedplant, plantZ)
            %     inputdata = myCosine(Test_time,Ts,Test_frequency,Test_amplitude);
            % elseif isequal(Designedplant, plantX)
            %     inputdata = mySine(Test_time,Ts,Test_frequency,Test_amplitude);
            % end
            sim('Block_fixedQ_Ccode_ideal');
            allrms(count_L,count_pole) = rms(ek);
            allmax(count_L,count_pole) = max(abs(ek));
            Test_inputdata = inputdata;
        end
    end
    
    csvwrite('C:/Users/dalek669528/Desktop/CNC_proj/MATLAB/data/' + date + '/' + Designedplant.name + '.csv', [reshape(poles, [length(poles), 1]), reshape(allQ, [size(allQ,1), size(allQ,3)])', allmax', allrms']);
    
    
    
    
    clear allFC fi k L Lq pole q Tcoe0_t Tcoe1_t Tcoe2_t Tphi tout
    % save("Z_chirp1hz_lengthtest.mat",'Designedplant','Q_length','poles','allfitness','allmax','allQ','allrms','allstatus','DesignedF','Reference_inputdata','Test_inputdata');
    % clc;
    % Simulation for Fixed parameters on Z axis
    % Expanded State Space with LFT
end

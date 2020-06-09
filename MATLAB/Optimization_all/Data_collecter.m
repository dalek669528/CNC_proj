clc;clear;close all;

date = "0609";
if ~exist("./data/" + date, 'dir')
   mkdir("./data/" + date)
   mkdir("./data/" + date + "/X")
   mkdir("./data/" + date + "/Z")
end
all_length = 15;
all_pole = 0.301:0.001:0.6;
% all_pole = 0.4952;
DesignT_change_times = 1;
randn_rate = 0.02;
% cvx_lims = [9e-2, 6e-2, 3e-2, 1e-2, 8e-3, 5e-3, 3e-3, 5e-4, 1e-4, 5e-5, 1e-5, 1e-6];
cvx_lims = 1e-5;

iter_times = length(all_length)*length(all_pole)*DesignT_change_times*length(cvx_lims);
checkpoint = iter_times;
max_err_X = zeros(1, iter_times);
max_err_Z = zeros(1, iter_times);
max_mid_err_X = zeros(1, iter_times);
max_mid_err_Z = zeros(1, iter_times);
rmse_err_X = zeros(1, iter_times);
rmse_err_Z = zeros(1, iter_times);
rmse_mid_err_X = zeros(1, iter_times);
rmse_mid_err_Z = zeros(1, iter_times);

DesignedF_Z_ori = [
    0.1     1000;
    1       100;
    10      10;
    100     1.3;
    347     0.831;
    509     0.4;
    700     0.11;
    1000    0.1;
    1400    0.1;
    1800    0.09;
    2150    0.09;
    2500    0.08;
    3000    0.083
    ];   

DesignedF_X_ori = [
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
ID_Model;
idx = 1;
for cvx_lim = cvx_lims
for Q_length = all_length
for poles = all_pole
for rand_times = 1:DesignT_change_times
    error = 0;
    for Designedplant = [plantZ, plantX]
        %% Design reference input
        Reference_time = 15; % (sec)
        Reference_amplitude = 50; % (mm)
        Reference_frequency = 2*pi*1/15;
        chirp_min = 0.01; % the lowest frequency
        chirp_max = 1.06; % the highest frequency
        inputdata = myChirpsine(Reference_time,Ts,chirp_min,chirp_max);
        Reference_inputdata = inputdata;

        %% Design the Frequency constraints
        if isequal(Designedplant, plantZ)
            DesignedF = DesignedF_Z_ori;
        elseif isequal(Designedplant, plantX)
            DesignedF = DesignedF_X_ori;
        end
        [Ncc,~] = size(DesignedF);
%         DesignedF(:,2) = abs(randn([Ncc, 1]).*DesignedF_ori(:,2)*randn_rate + DesignedF_ori(:,2));

        %% Solve the QCQP problem
        allfitness = zeros(Ncc,length(Q_length),length(poles));
        allQ = zeros(max(Q_length),length(Q_length),length(poles));
        allstatus = cell(length(Q_length),length(poles)); % cell array which can store any type of data
        allrms = zeros(2,length(Q_length),length(poles));
        allmax = zeros(2,length(Q_length),length(poles));
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
        %  use cvx to solve qcqp %
                tic;
                cvx_begin quiet
%                 args = [ cvx_lim, cvx_lim*5, cvx_lim*10];
%                 cvx_precision (args)
                cvx_precision low
                    variable q(L)
                    minimize((1.0/2.0) * quad_form(q,(2.0)*Tcoe2_t) + (-2*Tcoe1_t)*q + Tcoe0_t)
                    subject to
                    for fi = 1:allFC.Nc
                       (1.0/2.0) * quad_form(q,(2.0)*allFC.Omega(:,:,fi)) + allFC.beta(:,fi).'*q + allFC.gamma(fi) <= 0;
%                        (1.0/2.0) * 10^(-5) * quad_form(q,(2.0)*allFC.Omega(:,:,fi)) + 10^(-5) * allFC.beta(:,fi).'*q + 10^(-5) * allFC.gamma(fi) <= 0;
                    end
                cvx_end
                ggg = toc;

                fprintf("pole = %.3f, L = %d, status: %s, cvx: %f\n", pole, L, cvx_status, cvx_lim);

        % save data to txt %
                numberOfConstraints = allFC.Nc;
%                 save(".\data\numberOfConstraints_.txt", 'numberOfConstraints', '-ascii');
%                 save(".\data\Q_length_.txt", 'Q_length', '-ascii');
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
%                     save(".\data\P_matrix\P_matrix_" + fi + "_.txt", 'P_matrix', '-ascii');
%                     save(".\data\q_vector\q_vector_" + fi + "_.txt", 'q_vector', '-ascii');
%                     save(".\data\r_scalar\r_scalar_" + fi + "_.txt", 'r_scalar', '-ascii');
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
                    if isequal(Designedplant, plantZ)
                        inputdata = myCosine(Test_time,Ts,Test_frequency,Test_amplitude);
                    elseif isequal(Designedplant, plantX)
                        inputdata = mySine(Test_time,Ts,Test_frequency,Test_amplitude);
                    end
                    sim('Block_fixedQ_Ccode_ideal');
                    allrms(f,count_L,count_pole) = rms(ek);
                    allmax(f,count_L,count_pole) = max(abs(ek));
                end
                Test_inputdata = inputdata;
            end
        end
         clear allFC f fi k L Lq pole q Tcoe0_t Tcoe1_t Tcoe2_t Tphi tout
        % save("Z_chirp1hz_lengthtest.mat",'Designedplant','Q_length','poles','allfitness','allmax','allQ','allrms','allstatus','DesignedF','Reference_inputdata','Test_inputdata');
        % clc;
        % Simulation for Fixed parameters on Z axis
        % Expanded State Space with LFT
        if isequal(Designedplant, plantZ)
            dataZ = struct('name',"Z regular",...
                'pole',poles,...
                'Qnum',Qnum,...
                'rk',rk,...
                'yk',yk,...
                'ek',ek,...
                'uk',uk,...
                'DesignedF',DesignedF...
            );
        elseif isequal(Designedplant, plantX)
            dataX = struct('name',"X regular",...
                'pole',poles,...
                'Qnum',Qnum,...
                'rk',rk,...
                'yk',yk,...
                'ek',ek,...
                'uk',uk,...
                'DesignedF',DesignedF...          
            );
        end
        if ~(cellfun(@isequal, allstatus, {'Solved'}))
            disp('Not sloved')
            Test_time = 15*2;
            Test_amplitude = 50;
            Test_frequency = 2*pi*1/15;
            error = 1;
        end
    end
    
%     max_err_X = [];
%     max_err_Z = [];
%     max_mid_err_X = [];
%     max_mid_err_Z = [];
%     rmse_err_X = [];
%     rmse_err_Z = [];
%     rmse_mid_err_X = [];
%     rmse_mid_err_Z = [];
    if error == 0
        X_rk = dataX.rk;
        X_yk = dataX.yk;
        X_ek = dataX.ek;
        X_uk = dataX.uk;
        X_pole = dataX.pole;
        X_Qnum = dataX.Qnum;
        X_DesignedF = dataX.DesignedF;
        save("./data/" + date + "/X/X_" + string(sprintf('%.5d', idx)) + ".mat",'X_rk', 'X_yk', 'X_ek', 'X_uk', 'X_pole', 'X_Qnum', 'X_DesignedF')
        Z_rk = dataZ.rk;
        Z_yk = dataZ.yk;
        Z_ek = dataZ.ek;
        Z_uk = dataZ.uk;
        Z_pole = dataZ.pole;
        Z_Qnum = dataZ.Qnum;
        Z_DesignedF = dataZ.DesignedF;
        save("./data/" + date + "/Z/Z_" + string(sprintf('%.5d', idx)) + ".mat",'Z_rk', 'Z_yk', 'Z_ek', 'Z_uk', 'Z_pole', 'Z_Qnum', 'Z_DesignedF')
    end
    t = 0:Ts:Test_time;
    oneperiod = Test_time/4/Ts+1:Test_time/Ts-Test_time/4/Ts;
    max_err_Z(idx) = max(dataZ.ek)*1000;
    max_mid_err_Z(idx) = max(dataZ.ek(oneperiod))*1000;
    rmse_err_Z(idx) = rms(dataZ.ek)*1000;
    rmse_mid_err_Z(idx) = rms(dataZ.ek(oneperiod))*1000;
    max_err_X(idx) = max(dataX.ek)*1000;
    max_mid_err_X(idx) = max(dataX.ek(oneperiod))*1000;
    rmse_err_X(idx) = rms(dataX.ek)*1000;
    rmse_mid_err_X(idx) = rms(dataX.ek(oneperiod))*1000;
    idx = idx + 1;
end
end
end
end
figure('Name' ,'Checkpoint ' + string(idx), 'Position' ,[360 150 1200 780]);
subplot(2,2,1);hold on;grid on;
plot(max_mid_err_X,'r');
axis auto;
title('max(onewave) X ' + string(min(max_mid_err_X)));

subplot(2,2,3);hold on;grid on;
plot(rmse_mid_err_X,'r');
axis auto;
title('rmse(onewave) X ' + string(min(rmse_mid_err_X)));

subplot(2,2,2);hold on;grid on;
plot(max_mid_err_Z,'b');
axis auto;
title('max(onewave) Z ' + string(min(max_mid_err_Z)));

subplot(2,2,4);hold on;grid on;
plot(rmse_mid_err_Z,'b');
axis auto;
title('rmse(onewave) Z ' + string(min(rmse_mid_err_Z)));


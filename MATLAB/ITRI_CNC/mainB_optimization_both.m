clc;clear;close all;
spiral_data = load('./path/itrispiral.mat');
mycos_data = load('./path/itrimycos.mat');

%% Load the ID model from lib
ID_Model_ITRI;
Designedplant = plantZ_ITRI;
% referencedata = spiral_data.data.Z_spiral;%*Designedplant.count2round*Designedplant.round2mm;
testdata = mycos_data.data.Z_mycos;

%% Design reference input and Test input
% [Reference_inputdata, Reference_time]= CombineWithTime(referencedata,Ts);
Reference_time = 15;
chirp_min = 0.01;
chirp_max = 1.06;
% chirp_max = 0.56;
Reference_inputdata = myChirpsine(Reference_time,Ts,chirp_min,chirp_max);
Reference_inputdata(:,2) = Reference_inputdata(:,2)*Designedplant.mm2round*Designedplant.round2count;

[Test_inputdata, Test_time]= CombineWithTime(testdata,Ts); 
% Test_time = 15;
% Test_amplitude = 50;
% Test_frequency = 2*pi*2/15;
% Test_inputdata = myCosine(Test_time,Ts,Test_frequency,Test_amplitude);
% Test_inputdata(:,2) = Test_inputdata(:,2)*Designedplant.mm2round*Designedplant.round2count;

%% Design the Frequency constraints
DesignedF = [
            0.2     1000;
            2       100;
            20      10;
            100     1.1;
            200     0.8;
            300     0.7;
            400     0.6;
            500     0.5;
            600     0.4;
            700     0.3;
            800     0.2;
            900     0.1;
            1000    0.1;        
            1300    0.1;
            1500    0.1; 
            1800    0.1;
            2000    0.1;
            2500    0.1;
            3000    0.1            
            ];  
% DesignedF = [
%             0.1     1000;
%             1       100;
%             10      10;
%             100     1.1;
%             300     0.5;
%             400     0.33;
%             500     0.25;
%             600     0.2;
%             700     0.17;
%             800     0.14;
%             900     0.125;
%             1000    0.11;
%             1300    0.107;        
%             1500    0.105;
%             1800    0.103; 
%             2000    0.1;
%             2300    0.097;
%             2500    0.095;
%             2800    0.093;
%             3000    0.09            
%             ]; 
[Ncc,~] = size(DesignedF);

%% Determine the length and pole
% all_pole = (25:5:75)/100;
all_pole = 0.45;
all_length = 8:1:18;
% all_length = 12;
        
%% Solve the QCQP problem
allfitness = zeros(Ncc,length(all_length),length(all_pole));
allQ = zeros(max(all_length),length(all_length),length(all_pole));
allstatus = cell(length(all_length),length(all_pole));
allrms = zeros(2,length(all_length),length(all_pole));
allmax = zeros(2,length(all_length),length(all_pole));
count_pole = 0;

for pole = all_pole
    count_pole = count_pole + 1;
% Do Coprime Factorization  from lib
    plant = CoprimeFactorizationSS(Designedplant,pole);
    plant = LinearFractionalTransformation(plant);
    sim('Optimization_getVZ');
    
    count_L = 0;
    for L = all_length
        count_L = count_L + 1;
% Time-domaon coefficient %
        Tphi = zeros(L,1);
        Tcoe2_t = zeros(L,L);
        Tcoe1_t = zeros(1,L);
        Tcoe0_t = 0;
        for k = 1:length(vk)
            Tphi = [vk(k) ;Tphi(1:L-1)];
            Tcoe2_t = Tcoe2_t + Tphi*Tphi.';
            Tcoe1_t = Tcoe1_t + zk(k)*Tphi.';
            Tcoe0_t = Tcoe0_t + zk(k)^2;
        end
        Tcoe2_t = Tcoe2_t * 0.0001;
        Tcoe1_t = Tcoe1_t * 0.0001;
        Tcoe0_t = Tcoe0_t * 0.0001;
        
% Frequency domain coefficient %
        allFC = gConstraintFun(DesignedF,L,plant);
    
%  use cvx to solve qcqp %
        cvx_begin quiet
%         cvx_precision low
            variable q(L)
            minimize(quad_form(q,Tcoe2_t)-2*Tcoe1_t*q + Tcoe0_t)
            subject to
            for fi = 1:allFC.Nc
                quad_form(q,allFC.Omega(:,:,fi))+ allFC.beta(:,fi).'*q + allFC.gamma(fi) <= 0;
            end
        cvx_end
        fprintf("pole = %.3f, L = %d, status: %s\n", pole, L,cvx_status);

% record data %
        allfitness(:,count_L,count_pole) = CheckFitnessFun(q,allFC);
        allstatus(count_L,count_pole) = {cvx_status};
        allQ(1:L,count_L,count_pole) = q;

% perfomance test %
        [ccA,ccB,ccC,ccD] = LFTExpandedSS(plant.Jy,q);
        C = ss(ccA,ccB,ccC,ccD,Ts);
        Qnum = q;
        Lq = L; % for simulink
        inputdata = Test_inputdata;
        sim('ITRI_Block_fixedQ_Ccode');
        allrms(1,count_L,count_pole) = rms(ek_exp);%*plant.mm2round*plant.round2count;
        allrms(2,count_L,count_pole) = rms(ek_exp)*plant.count2round*plant.round2mm;
        allmax(1,count_L,count_pole) = max(abs(ek_exp));%*plant.mm2round*plant.round2count;
        allmax(2,count_L,count_pole) = max(abs(ek_exp))*plant.count2round*plant.round2mm;
%         PlotMagnitudeWithConstraints('',C*plant.v2p,DesignedF);
    end
end
clear allFC f fi k L pole q Tcoe0_t Tcoe1_t Tcoe2_t Tphi tout %Lq inputdata

%%
% save('X_itrichirp_regular_v1.mat','Designedplant','all_length','all_pole','allfitness','allmax','allQ','allrms','allstatus','DesignedF','Reference_inputdata','Test_inputdata');
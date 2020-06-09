% Simulation for Fixed parameters on Z axis
% Expanded State Space with LFT
clc;clear;close all;
mycosdata = load('./path/itrimycos.mat');
load('./path/path3030.mat');

%% Load the ID model from lib
ID_Model_ITRI;
pole = 0.45;
plant = plantZ_ITRI;

%% Do Coprime Factorization  from lib
plant = CoprimeFactorizationSS(plant,pole);

%% Do Coprime Factorization via LFT from lib
plant = LinearFractionalTransformation(plant);

%% Determine the FIR parameters
% regular with chirp
% Qnum = [520.834317367458,72.3760875691095,206.540186517313,136.302013214289,339.420144351987,299.799326355582,314.955904831065,279.207588415117,321.964942320754,261.979588287706,180.389773201754,-45.1045838142933,-207.809678492573,-162.447211518657,-770.056299825798];
Qnum = [884.956187600194,-237.650887872230,25.2002633846358,225.998112609307,-108.428184518527,69.8382065610213,76.5980590386770,-99.7888652893205,176.863086182581,166.811578384177,-30.4802889138195,594.924959164552];
Lq = length(Qnum);

%% Check the constraints fitness
DesignedF = [
            0.1     1000;
            1       100;
            10      10;
            100     1.1;
            300     0.8;
            400     0.7;
            500     0.6;
            600     0.5;
            700     0.4;
            800     0.3;
            900     0.2;
            1000    0.1;
            1300    0.1;        
            1500    0.1;
            1800    0.1; 
            2000    0.1;
            2500    0.1;
            3000    0.1            
            ];   

allFrequencyconstraints = gConstraintFun(DesignedF,Lq,plant);
allfitness = CheckFitnessFun(Qnum.',allFrequencyconstraints);

%% Decide the Controller
[ccA,ccB,ccC,ccD] = LFTExpandedSS(plant.Jy,Qnum);
CC = ss(ccA,ccB,ccC,ccD,Ts);

OLoop = minreal(CC*plant.v2p);
SLoop = minreal(1/minreal(1+OLoop));
CLoop = minreal(OLoop*SLoop);

%% Run Simulation on Simulink
% Test_time = 30;
% Test_amplitude = 50;
% Test_frequency = 2*pi*1/15;
% inputdata = myDesignedCosine(Test_time,Ts,Test_frequency,Test_amplitude);
% inputdata = myCosine(Test_time,Ts,Test_frequency,Test_amplitude);
% [inputdata, Test_time]= CombineWithTime(mycosdata.data.X_mycos*plant.count2round*plant.round2mm,Ts);
% [inputdata, Test_time]= CombineWithTime(mycosdata.data.X_mycos,Ts);
[inputdata, Test_time]= CombineWithTime(pathdata.Zin*plant.mm2round*plant.round2count,Ts);

sim('ITRI_Block_fixedQ_Ccode');


%% Show Frequency-domain Results
% close all;
% PlotMagnitudeWithConstraints('|L|',OLoop,[]);
% PlotMagnitudeWithConstraints('|L|',OLoop,DesignedF);
% PlotMagnitudeWithConstraints('|1/(1+L)|',SLoop,[]);
% PlotMagnitudeWithConstraints('|L/(1+L)|',CLoop,[]);

% margin(OLoop);grid on;
% margin(SLoop);grid on;
% margin(CLoop);grid on;

%% Show Time-domain Results
t = 0:Ts:Test_time;
oneperiod = round(Test_time/4/Ts+1):round(Test_time/Ts-Test_time/4/Ts);

PlotTimeDomainResult_ITRI('Fixed parameters - Z - expanded C code',rk_exp,yk_exp,ek_exp,t,plant);

max_rmse_exp = max(ek_exp);
rmse_exp = rms(ek_exp);
rmse_mid_exp = rms(ek_exp(oneperiod));
fprintf("(expanded C code)\nL: %d\nMax error: %.6f counts\nRMS(all): %.6f counts\nRMS(onewave): %.6f counts\n", Lq, max_rmse_exp, rmse_exp, rmse_mid_exp);

%%
% myf = fopen('errorZ.txt','w');
% for i = 1:length(ek_exp)
%     fprintf(myf,"%d\r\n",round(ek_exp(i)));
% end
% fclose(myf);
%%
% myf = fopen('ukZ.txt','r');
% Cuk = fscanf(myf,"%lf");
% fclose(myf);
% dif = uk_exp-Cuk;
% fclose('all');
% Simulation for PD+FF of ITRI
clc;clear;close all;
load('./path/reference3775.mat');
load('./path/reference0.mat');
load('./path/path3030.mat');
mycosdata = load('./path/itrimycos.mat');

% Load the ID model from lib
ID_Model_ITRI;
plant = plantZ_ITRI.v2p;

% rising time: 0.04 -> robus:0.7 -> badwidth:70 -> phase margin:70
% Design the X controller parameter
Kp_x = 65.4554;
Kd_x = 0.090865;
F_x = 1;
Cx = Kp_x + Kd_x * 1/Ts *(z-1)/z;
FFx = F_x*(z-1)/Ts/z;

% Design the Z controller parameter
Kp_z = 66.0566;
Kd_z = 0.14788;
F_z = 1;
Cz = Kp_z + Kd_z * 1/Ts *(z-1)/z;
FFz = F_z*(z-1)/Ts/z;

% Run simulink simulation
% Reference_time = (length(reference3775.X_cos)*2-1)*Ts;
% inputdata_x = [(0:Ts:Reference_time).' [reference3775.X_cos;reference3775.X_cos]];
% inputdata_z = [(0:Ts:Reference_time).' [reference3775.Z_cos;reference3775.Z_cos]];
% [inputdata_x, Reference_time]= CombineWithTime(mycosdata.data.X_mycos,Ts);
% [inputdata_z, ~]= CombineWithTime(mycosdata.data.Z_mycos,Ts);
[inputdata_x, Reference_time]= CombineWithTime(pathdata.Xin*plantX_ITRI.mm2round*plantX_ITRI.round2count,Ts);
[inputdata_z, ~]= CombineWithTime(pathdata.Zin*plantZ_ITRI.mm2round*plantZ_ITRI.round2count,Ts);
% Reference_time = 15;
% Reference_amplitude = 50;
% Reference_frequency = 2*pi*2/15;
% inputdata = myCosine(Reference_time,Ts,Reference_frequency,Reference_amplitude);
% inputdata_x(:,1) = inputdata(:,1);
% inputdata_x(:,2) = round(inputdata(:,2)*plantX_ITRI.mm2round*plantX_ITRI.round2count);
% inputdata_z(:,1) = inputdata(:,1);
% inputdata_z(:,2) = round(inputdata(:,2)*plantZ_ITRI.mm2round*plantZ_ITRI.round2count);
sim('ITRI_Block_PDFF_Ccode_XZ');

%% Plot bode diagram
% close all;
% figure;bode((FFx+Cx)*plantX.v2p/(1+Cx*plantX.v2p),0.01:0.1:pi/Ts);grid on;
% figure;bode((FFz+Cz)*plantZ_ITRI.v2p/(1+Cz*plantZ_ITRI.v2p),0.01:0.1:pi/Ts);grid on;
% figure;bode((FFz+Cz)*plantZ_ITRI.v2p,0.01:0.1:pi/Ts);grid on;
% figure;bode((FFx+Cx)*plantX_ITRI.v2p,0.01:0.1:pi/Ts);grid on;
% figure;margin((FFz+Cz)*plantZ.v2p/(1+Cz*plantZ.v2p));grid on;

%% Plot signal messgaes
t = 0:Ts:Reference_time;

PlotTimeDomainResult_ITRI('Z-axis',rk_Z,yk_Z,ek_Z,t,plantZ_ITRI);
max_rmse_Z = max(ek_Z);
RMSE_Z_count = rms(ek_Z);
RMSE_Z_mm = RMSE_Z_count*plantZ_ITRI.count2round*plantZ_ITRI.round2mm;
fprintf("Z:%f\nMax error: %f\nRMS(all): %f (%f mm)\n",F_z,max_rmse_Z, RMSE_Z_count,RMSE_Z_mm);

PlotTimeDomainResult_ITRI('X-axis',rk_X,yk_X,ek_X,t,plantX_ITRI);
max_rmse_X = max(ek_X);
RMSE_X_count = rms(ek_X);
RMSE_X_mm = RMSE_X_count*plantX_ITRI.count2round*plantX_ITRI.round2mm;
fprintf("X:%f\nMax error: %f\nRMS(all): %f (%f mm)\n",F_x,max_rmse_X, RMSE_X_count,RMSE_X_mm);
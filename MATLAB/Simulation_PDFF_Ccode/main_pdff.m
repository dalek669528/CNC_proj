% Simulation for PD+FF on XZ plant
clc;clear;close all;

% Load the ID model from lib
ID_Model;

plant  = plantZ.v2p;
% rising time: 0.04 -> robus:0.7 -> badwidth:70 -> phase margin:70
% Design the X controller parameter
Kp_x = 460;
Kd_x = 4;
F_x = 5.975;
Cx = Kp_x + Kd_x * 1/Ts *(z-1)/z;
FFx = F_x*(z-1)/Ts/z;

% Design the Z controller parameter
Kp_z = 355;
Kd_z = 4;
F_z = 4.955;
Cz = Kp_z + Kd_z * 1/Ts *(z-1)/z;
FFz = F_z*(z-1)/Ts/z;

% Run simulink simulation
Reference_time = 15*2;
Reference_amplitude = 50;
Reference_frequency = 2*pi*1/15;

inputdata_x = myDesignedCosine(Reference_time,Ts,Reference_frequency,Reference_amplitude);
% inputdata_x = mySine(Reference_time,Ts,Reference_frequency,Reference_amplitude);
% inputdata_x = myCosine(Reference_time,Ts,Reference_frequency,Reference_amplitude);
inputdata_z = inputdata_x;
sim('Block_PDFF_XZ_Ccode')

%% Plot bode diagram
% close all;
% figure;bode((FFx+Cx)*plantX.v2p/(1+Cx*plantX.v2p),0.01:0.1:pi/Ts);grid on;
% figure;bode((FFz+Cz)*plantZ.v2p/(1+Cz*plantZ.v2p),0.01:0.1:pi/Ts);grid on;
% figure;bode((FFz+Cz)*plantZ.v2p,0.01:0.1:pi/Ts);grid on;
% figure;margin((FFz+Cz)*plantZ.v2p/(1+Cz*plantZ.v2p));grid on;

%% Plot signal messgaes
t = 0:Ts:Reference_time;
oneperiod = Reference_time/4/Ts+1:Reference_time/Ts-Reference_time/4/Ts;

PlotTimeDomainResult('Z-axis',rk_Z,yk_Z,ek_Z,t);
max_rmse_Z = max(ek_Z)*1000;
RMSE_Z = rms(ek_Z)*1000;
RMSE_mid_Z = rms(ek_Z(oneperiod))*1000;
fprintf("Z:\nMax error: %f\nRMS(all): %f\nRMS(onewave): %f\n",max_rmse_Z, RMSE_Z, RMSE_mid_Z);

PlotTimeDomainResult('X-axis',rk_X,yk_X,ek_X,t);
max_rmse_X = max(ek_X)*1000;
RMSE_X = rms(ek_X)*1000;
RMSE_mid_X = rms(ek_X(oneperiod))*1000;
fprintf("X:\nMax error: %f\nRMS(all): %f\nRMS(onewave): %f\n",max_rmse_X, RMSE_X, RMSE_mid_X);

%%
% data = struct('name',"X pdff",'rk',rk_X,'yk',yk_X,'ek',ek_X,'FF',F_x,'Kp',Kp_x,'Kd',Kd_x);
% save("X_pdff_1.mat",'data')
% data = struct('name',"Z pdff",'rk',rk_Z,'yk',yk_Z,'ek',ek_Z,'FF',F_z,'Kp',Kp_z,'Kd',Kd_z);
% save("Z_pdff_1.mat",'data')

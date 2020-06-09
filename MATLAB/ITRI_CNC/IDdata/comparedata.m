clc; clear; close all;
G3C = load('G3C.mat');
G2C = load('G2C.mat');
G3S = load('G3S.mat');
G2S = load('G2S.mat');

Xin = G2C.mmdata.X_rk;%*plantX_ITRI.count2round*plantX_ITRI.round2mm;
Zin = G2C.mmdata.Z_rk;%*plantZ_ITRI.count2round*plantZ_ITRI.round2mm;
t = (0:length(Xin)-1)*G3C.mmdata.Ts;

% range = 100:6774; radius = 30; center = [30 0]; feedrate = 2000;% G3C
range = 100:6201; radius = 30; center = [30 0]; feedrate = 1300;% G2C

% pathdata = struct('Xin',Xin,'Zin',Zin);save('path3030.mat','pathdata');
% plot(-Zin,Xin);axis square;title('Reference path (G03) (CCW)');xlabel('(mm)');ylabel('(mm)')
% subplot(211);plot(t,Xin);axis([-inf inf -30 30]);title('X path');ylabel('(mm)');
% subplot(212);plot(t,Zin);axis([-inf inf 0 -60]);title('Z path');xlabel('time(s)');ylabel('(mm)');


[G3C_M,G3C_m,G3C_Roundness] = PlotRealCircle('G3C ITRI',G3C.mmdata,range,radius,center);
% [G2C_M,G2C_m,G2C_Roundness] = PlotRealCircle('G2C ITRI',G3C.mmdata,range,radius,center);

% [nctu_rmseX, nctu_maxeX, nctu_rmseZ, nctu_maxeZ] = PlotTwoError('nctu ITRI',nctuXZ.mmdata,range);


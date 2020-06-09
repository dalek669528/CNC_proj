% sampling time
Ts = 0.001;
z = tf('z',Ts);
Vmax = 300;  %unknown
%---------------------------- initial -----------------------------%
% % transfer function (x axis) (cps 2 cps)
% Pv_x = tf([0,0.184613578780250,-0.357467968465759,0.173016592756523,0],...
%           [1,-2.81134086427335,2.66789703375917,-0.889675380505102,0.0332810815356037],z,'Inputdelay',2);
%     
% % transfer function (z axis) (cps 2 cps)
% Pv_z = tf([0,0.173227652422342,-0.332122085591785,0.159061916413243,0],...
%           [1,-2.60556002424451,2.06540929919667,-0.299143719032096,-0.160538316460838],z,'Inputdelay',2);
%---------------------------- initial -----------------------------%

% transfer function (x axis) (cps 2 cps)
Pv_x = tf([0,0.0193,0.3593,-0.3301],...
          [1,-1.9914,1.6077,-0.7181,0.1505],z,'Inputdelay',0);
    
% transfer function (z axis) (cps 2 cps)
Pv_z = tf([0,0.1176,-0.0697,-0.0097],...
          [1,-2.3436,1.9912,-0.6962,0.0866],z,'Inputdelay',0);

integrator = tf([Ts,0],[1,-1],z);

X_count2round = 1/8388608;
Z_count2round = 1/8388608;
X_mm2round = 1/10;
Z_mm2round = 1/12;

% transfer function (z axis) (cps 2 count)
Pz = Pv_z*integrator;
% transfer function (x axis) (cps 2 count)
Px = Pv_x*integrator;

%%
plantZ_ITRI = struct('v2p',Pz,'v2v',Pv_z,'Ts',Ts,'count2round',Z_count2round,'mm2round',Z_mm2round,'round2count',1/Z_count2round,'round2mm',1/Z_mm2round);
plantX_ITRI = struct('v2p',Px,'v2v',Pv_x,'Ts',Ts,'count2round',X_count2round,'mm2round',X_mm2round,'round2count',1/X_count2round,'round2mm',1/X_mm2round);
clear Pv_x Pv_z Px Pz integrator X_count2round Z_count2round X_mm2round Z_mm2round;

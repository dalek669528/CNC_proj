% sampling time
Ts = 0.001;
z = tf('z',Ts);
Vmax = 300;

% transfer function (x axis) (rpm 2 rpm)
Pv_x = tf([0,0.032229639522932,0.060985247548587,-0.189017468410220,0.096991370106707],...
        [1,-3.084824116823246,3.534522236179820,-1.782062555742226,0.333547393790890],z);
    
% transfer function (z axis) (rpm 2 rpm)
Pv_z = tf([0,0.0973217504156460,-0.209580502231482,0.151454774135651,-0.0370139184925609],...
        [1,-3.44700230387150,4.55138855100522,-2.73113485001183,0.628907804339994],z);
    
integrator = tf([Ts,0],[1,-1],z);

rpm2mms_z = 12/60;
rpm2mms_x = 10/60;

% transfer function (z axis) (rpm 2 mms)
Pz = Pv_z*rpm2mms_z*integrator;
% transfer function (x axis) (rpm 2 mms)
Px = Pv_x*rpm2mms_x*integrator;

% uncertainty = 1; % no uncertainty

%%
plantZ = struct('v2p',Pz,'v2v',Pv_z,'Ts',Ts);
plantX = struct('v2p',Px,'v2v',Pv_x,'Ts',Ts);
clear Pv_x Pv_z Px Pz integrator uncertainty;

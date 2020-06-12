% Simulation for Fixed parameters on X axis
% Expanded State Space with LFT
clc;clear;close all;

%% Determine the FIR parameters
% Q is order N (length = L = N+1), 
% FIR Q = q0 + q1 z^-1 +... + qn z^-N
%       = (q0 z^N + q1 z^(N-1) + ... + qn)/z^N

% time only / L = 5
% Qnum = [22280.6409756861;12218.5694012669;2305.56982548507;-7046.68934227489;-15597.5733636341].';

% regular designed
% pole = 0.25;Qnum = [7720.25072178680;-1877.89641705759;3771.54152461758;1399.47065336347;3197.19987540252].';
% pole = 0.25;Qnum = [9089.40006386962;-2910.11887235360;3057.24480707413;2382.43724236762;1539.99205145195;1201.94768438137;1000.86899781272;-672.394450207280;383.824294416111;-864.435832579171].';
% pole = 0.25;Qnum = [7322.63635475436;-1433.04411453800;3112.85910050061;2034.32758486280;2443.07848612195;2980.08111169159;2468.58138976403;2398.39242960770;1274.99717639972;851.896688278147;-633.157559928313;-314.781808035084;-2964.76681553323;-475.013491670385;-4882.39911642268].';

% one peak
% pole = 0.25;Qnum = [7627.80160667853;-1917.64123260947;3752.89004427902;1373.24808969987;3246.89201911388].';
% pole = 0.25;Qnum = [8672.94312132963;-3206.59632277202;3813.96583553142;2735.63742473263;1464.88160296393;842.958629363714;408.646100614129;-692.172324507277;579.201438842517;-409.776721776658].';
% pole = 0.25;Qnum = [6773.99474567456;-2653.29466719424;3193.34755491539;2716.83457343600;2661.93821862870;2953.65755849207;2283.15480720508;2434.12455009330;1445.39194603424;986.567953251885;-884.265455196822;-874.109011603139;-3451.03450992730;337.542694864155;-3737.02501822945].';

% two
% pole = 0.25;Qnum = [7627.42963188591;-1917.01595573582;3752.90783217394;1372.60207496746;3247.27203973194].';
% pole = 0.25;Qnum = [7631.71982577424;-1927.94306704964;3691.31912121638;1423.25483399251;3264.44993324413;152.717696271927;103.175570217596;-28.8473052680319;68.1914816559489;-168.419766528264].';
% pole = 0.25;Qnum = [6363.48848233640;-2109.82853811019;3158.05324574143;2105.01178710208;3398.53893208179;2695.18036490339;2150.07872032454;2341.14751220413;1523.57332529841;1010.41679963907;-1543.14853929645;-106.016541223867;-3563.07542104437;5.98938650539695;-3241.81862652875].';

% % % new
% regular
% pole = 0.4352; Qnum = [8885.43106204057;-5053.40645444674;4563.93354254597;731.929927986032;586.133097341372].';
% pole = 0.4352; Qnum = [7725.47436870418;-2800.76936262095;3548.91562110505;1455.85047142394;1551.92560685381;2259.01570131598;1563.23963389991;1531.89199295571;431.816432445653;490.165544174720;-596.671338708589;173.846343973470;-3658.70856106041;1015.91649175504;-4995.63495385305].';
% notch
% pole = 0.4352; Qnum = [8723.45113192471;-4898.56030515831;4553.34552488401;557.177856361654;772.601425788690].';
% pole = 0.4352; Qnum = [7473.04093443492;-3560.91689973965;3344.74407264701;1512.08207161222;1405.87592271814;2116.43279429573;1351.87854986277;1379.05991114751;314.771935700433;497.292290224317;-533.575868853380;345.335872655316;-3440.74107426101;1735.76314107882;-4242.63151155773].';
% pole = 0.4512; Qnum = [8245.1631765846, -4335.4780233633, 3193.8850472014,  2472.2664581544,  1180.2966100668, 1916.7070906676,  1307.2304888973,  1207.8277786821, -756.8689804559,   986.7789200112, -1757.8461805759,  -87.7924366883, -2056.0527675153,   209.634839905 ,-2378.7260136405 ];
pole = 0.4512; Qnum = [8654.903  , -6973.163  ,   726.6929 ,  1748.4072 , -1208.1191 ,  1739.7069 ,  -103.84458, -2114.6292 ,  -553.1921 , 1104.388  , -3578.4001 ,  2252.313  , -3883.007  ,   672.8919 , -2579.314 ];
Lq = length(Qnum);

%% Load the ID model from lib
ID_Model;
plant = plantX;

%% Do Coprime Factorization  from lib
plant = CoprimeFactorizationSS(plant,pole);

%% Do Coprime Factorization via LFT from lib
plant = LinearFractionalTransformation(plant);

%% Check the constraints fitness
% new design (0.4352/5) X
DesignedF = [
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
allFrequencyconstraints = gConstraintFun(DesignedF,Lq,plant);
allfitness = CheckFitnessFun(Qnum.',allFrequencyconstraints);

% use this command to save constraints for LabVIEW
% Saveconstraints('X_Notch_15',allFrequencyconstraints);

%% Decide the Controller
[ccA,ccB,ccC,ccD] = LFTExpandedSS(plant.Jy,Qnum);
CC = ss(ccA,ccB,ccC,ccD,Ts);

OLoop = minreal(CC*plant.v2p);
SLoop = minreal(1/minreal(1+OLoop));
CLoop = minreal(OLoop*SLoop);

%% Run Simulation on Simulink
Test_time = 15*2;
Test_amplitude = 50;
Test_frequency = 2*pi*1/15;
% inputdata = myDesignedCosine(Test_time,Ts,Test_frequency,Test_amplitude);
inputdata = myCosine(Test_time,Ts,Test_frequency,Test_amplitude);
% inputdata = mySine(Test_time,Ts,Test_frequency,Test_amplitude);

% sim('Block_fixedQ_Ccode_expanded');

sim('Block_fixedQ_Ccode_ideal');
%% Show Frequency-domain Results
% close all;
% PlotMagnitudeWithConstraints('|L|',OLoop,[]);
PlotMagnitudeWithConstraints('|L|',OLoop,DesignedF);
% PlotMagnitudeWithConstraints('|1/(1+L)|',SLoop,[]);
% PlotMagnitudeWithConstraints('|L/(1+L)|',CLoop,[]);

% margin(OLoop);grid on;
% margin(SLoop);grid on;
% margin(CLoop);grid on;

%% Show Time-domain Results
t = 0:Ts:Test_time;
oneperiod = Test_time/4/Ts+1:Test_time/Ts-Test_time/4/Ts;

PlotTimeDomainResult('Fixed parameters - X',rk,yk,ek,t);

max_rmse_ss = max(ek)*1000;
rmse_ss = rms(ek)*1000;
rmse_mid_ss = rms(ek(oneperiod))*1000;
fprintf("(stste space)\nL: %d\nMax error: %.6f um\nRMS(all): %.6f um\nRMS(onewave): %.6f um\n", Lq, max_rmse_ss, rmse_ss, rmse_mid_ss);

%%  save data
% data = struct('name',"X regular",...
%               'pole',pole,...
%               'Qnum',Qnum,...
%               'rk',rk_exp,...
%               'yk',yk_exp,...
%               'ek',ek_exp,...
%               'DesignedF',DesignedF...          
%               );
% save("X_n_5_1.mat",'data')
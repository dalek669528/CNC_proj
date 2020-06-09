% Simulation for Fixed parameters on Z axis
% Expanded State Space with LFT
clc;clear;close all;

%% Determine the FIR parameters
% Q is order N (length = L = N+1), 
% FIR Q = q0 + q1 z^-1 +... + qn z^-N
%       = (q0 z^N + q1 z^(N-1) + ... + qn)/z^N

% time only
% pole = 0.5;Qnum = [10953.235095878084, 6011.039747600123, 1214.787993583139, -3282.0109200485404, -7387.394950718186];

% regular with chirp
% pole = 0.5;Qnum = [3973.87865353511;-2027.97337669182;2510.79198949738;-2436.73382336036;5489.35148861471].';
% pole = 0.5;Qnum = [8972.07278169804;-4921.16690120231;2686.33608940740;-224.009837396147;-274.740815006160;-338.013561236905;-355.073427910869;389.737625914920;537.012978367426;1045.17137711001].';
% pole = 0.5;Qnum = [4332.68720809435;-762.824007058915;3079.39627759621;-548.742069565578;583.662600399268;-139.519105490218;531.489726955241;465.655333713069;1073.86410215099;808.069244251544;1769.88829799787;469.778970291169;-433.911680836792;3754.02296145848;-7477.21403664050].';

% one notch with chirp
% pole = 0.5;Qnum = [8030.25776202586;-5734.92057533428;1602.41483518097;1306.54882991698;2304.24547745233].';
% pole = 0.5;Qnum = [9216.76432891263;-6099.97520285813;-1982.29321840382;4026.71260476584;5342.81962918442;-2167.70266710318;-3228.79685685113;1373.51948642062;3575.82988687572;-2541.48983947885].';
% pole = 0.5;Qnum = [15494.4463944286;-14114.9213344333;-7920.37792365175;16511.2553581611;10775.5649389209;-11425.6043387092;-13476.7694487267;11342.4966820056;14801.7475283938;-5858.10916478414;-15784.2783851814;7128.27375689334;7468.97090386769;-4692.71501710725;-2743.42941104026].';

% two notch with chirp
% pole = 0.5;Qnum = [4910.26376707223;84.4680515332290;-1711.63061502950;325.777727548168;3911.87690686493].';
% pole = 0.5;Qnum = [4938.36557085589;121.534491172765;-1766.31619897105;272.488748431873;3990.08409352744;43.3035716063373;-41.4081038884260;-52.5444224785757;31.6198981969520;-21.0322430837634].';
% pole = 0.5;Qnum = [4636.90215899127;-614.974954914681;-822.580116171607;1173.87010625088;2618.67126989000;256.809213517662;843.262738306668;489.027727900442;-47.9627099343472;1209.09267655744;-1723.31150164160;323.593944990977;633.504041229690;-41.1118049461247;-1423.80024796071].';

% % % new
% regular
% --> pole = 0.4952;Qnum = [3743.42046414929;-1843.07043055318;2492.03430308000;-2520.90147252994;5731.14592620173].';

% pole = 0.4952;Qnum = [4359.83206978553;-743.950262090048;3063.60418151860;-538.039851257320;575.319927204521;-123.379365953110;534.228637206103;460.216323213538;1072.93539216552;826.938754146939;1763.14067495987;497.782650416741;-405.489735135970;3689.90656854944;-7430.62801188454].';
% notch
% pole = 0.4952;Qnum = [3848.76700192144;-1948.17891885533;2521.74310409079;-2445.85314220169;5578.57064207358].';
%pole = 0.4952;Qnum = [4670.29807752366;-2915.06315138601;2492.97016243023;-736.557047422892;1800.56681396059;1362.26065656813;510.851838799402;966.655615129803;1601.77515826129;366.553977436919;-472.515804825815;1471.75525011512;-2359.70115655455;2396.48621981079;-3552.30879248213].';% C++ solved
pole = 0.4952;Qnum = [4670.21;-2916.04;2493.59;-737.039;1800.59;1362.52;511.299;967.203;1602.43;366.933;-471.903;1470.56;-2359.58;2397.5;-3554.25].';
Lq = length(Qnum);

%% Load the ID model from lib
ID_Model;
plant = plantZ;

%% Do Coprime Factorization  from lib
plant = CoprimeFactorizationSS(plant,pole);

%% Do Coprime Factorization via LFT from lib
plant = LinearFractionalTransformation(plant);

%% Check the constraints fitness
% new design (0.4952/5) Z       
DesignedF = [
            0.1     1000;
            1       100;
            10      10;
            100     1.1;
            300     0.8;
            400     0.56;
            500     0.39;
            600     0.26;
            700     0.16;
            800     0.09;
            900     0.07;
            1000    0.08;
            1300    0.1;        
            1500    0.09;
            1800    0.06; 
            2000    0.03;
            2500    0.07;
            3000    0.1            
            ]; 
allFrequencyconstraints = gConstraintFun(DesignedF,Lq,plant);
allfitness = CheckFitnessFun(Qnum.',allFrequencyconstraints);

% use this command to save constraints for LabVIEW
% Saveconstraints('Z_Notch_15',allFrequencyconstraints);

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

sim('Block_fixedQ_Ccode_expanded');

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
oneperiod = Test_time/4/Ts+1:Test_time/Ts-Test_time/4/Ts;

PlotTimeDomainResult('Fixed parameters - Z - state space',rk,yk,ek,t);

max_rmse_ss = max(ek)*1000;
rmse_ss = rms(ek)*1000;
rmse_mid_ss = rms(ek(oneperiod))*1000;
fprintf("(stste space)\nL: %d\nMax error: %.6f um\nRMS(all): %.6f um\nRMS(onewave): %.6f um\n", Lq, max_rmse_ss, rmse_ss, rmse_mid_ss);

%%  save data
% data = struct('name',"Z regular",...
%               'pole',pole,...
%               'Qnum',Qnum,...
%               'rk',rk_exp,...
%               'yk',yk_exp,...
%               'ek',ek_exp,...
%               'DesignedF',DesignedF...
%               );
% save("Z_n_5_1.mat",'data')          
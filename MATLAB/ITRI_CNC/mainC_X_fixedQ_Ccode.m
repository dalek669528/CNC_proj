% Simulation for Fixed parameters on X axis
% Expanded State Space with LFT
clc;clear;close all;
mycosdata = load('./path/itrimycos.mat');
load('./path/path3030.mat');

%% Load the ID model from lib
ID_Model_ITRI;
pole = 0.5;
plant = plantX_ITRI;

%% Do Coprime Factorization  from lib
plant = CoprimeFactorizationSS(plant,pole);

%% Do Coprime Factorization via LFT from lib
plant = LinearFractionalTransformation(plant);

%% Determine the FIR parameters
% regular with chirp
% Qnum = [459.894575946634,-16.7276801647275,237.320750790429,-4.60207528387945,199.839040629263,174.674836860375,201.223687077333,139.281611653196,191.204486586651,159.528088250819,175.512613312057,-9.56756427385017,-159.586683824721,94.6386097953774,-764.850435983925];
% Qnum = [932.633538932108;-505.043820594633;320.993766069384;-13.6823147892448;1.56859894131693;78.7364827789356;-147.954840493425;28.6806181333443;-46.8779431343784;89.0077187114210;-156.090152279880;311.679743829099;56.0491034530710;-99.3050331478834;225.314491241860].';
% Qnum = [455.217632699638;-275.845978412654;16.4952079939694;-91.5944146982562;6.54836094421551;6.25366942189307;85.5123441602310;112.511580362748;269.813471021598;-48.9744028124197;529.607264077809].';
Qnum = [797.514331990504,-485.232610846030,262.856414138412,24.3359377967718,-149.423463056821,175.247356782652,-166.053327628022,61.0772708379591,32.3954543787851,126.465703010449,-51.3265230853852,444.837987690054];
% Qnum = 300*ones(1,15);
Lq = length(Qnum);

%% Check the constraints fitness
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
% inputdata = mySine(Test_time,Ts,Test_frequency,Test_amplitude);
% inputdata(:,2) = inputdata(:,2)*plant.mm2round*plant.round2count;
% [inputdata, Test_time]= CombineWithTime(round(mycosdata.data.X_mycos),Ts);
[inputdata, Test_time]= CombineWithTime(pathdata.Xin*plant.mm2round*plant.round2count,Ts);
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

PlotTimeDomainResult_ITRI('Fixed parameters - X - expanded C code',rk_exp,yk_exp,ek_exp,t,plant);

max_rmse_exp = max(ek_exp);
rmse_exp = rms(ek_exp);
rmse_mid_exp = rms(ek_exp(oneperiod));
fprintf("(expanded C code)\nL: %d\nMax error: %.6f counts\nRMS(all): %.6f counts\nRMS(onewave): %.6f counts\n", Lq, max_rmse_exp, rmse_exp, rmse_mid_exp);

%%
% myf = fopen('errorX.txt','w');
% for i = 1:length(ek_exp)
%     fprintf(myf,"%d\r\n",round(ek_exp(i)));
% end
% fclose(myf);
%%
% myf = fopen('ukX.txt','r');
% Cuk = fscanf(myf,"%lf");
% fclose(myf);
% dif = uk_exp-Cuk;
% fclose('all');

% Simulation for Fixed parameters on X axis
% Expanded State Space with LFT
clc;clear;close all;
client = tcpclient('localhost', 50000);
csv_table = csvread('C:\Users\dalek669528\Desktop\CNC_proj\MATLAB\data\Generate\temp.csv');
data_amount = length(csv_table);
% data_amount = 10;
% max_ss_all = zeros(data_amount,1);
max_ss_one = zeros(data_amount,1);
% rmse_ss_all = zeros(data_amount,1);
rmse_ss_one = zeros(data_amount,1);

for idx = 1:data_amount
    pole = csv_table(idx, 1);
    Qnum = csv_table(idx, 2:16);
%     DesignedF = [0.1, 1, 10, 100, 300, 400, 500, 600, 700, 800, 900, 1000, 1300, 1500, 1800, 2000, 2500, 3000];
%     DesignedF = [DesignedF; csv_table(idx, 3:20)]';
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
    %% Determine the FIR parameters
    % Q is order N (length = L = N+1), 
    % FIR Q = q0 + q1 z^-1 +... + qn z^-N
    %       = (q0 z^N + q1 z^(N-1) + ... + qn)/z^N
%     pole = 0.4352; 
%     QnuRm = [7473.04093443492;-3560.91689973965;3344.74407264701;1512.08207161222;1405.87592271814;2116.43279429573;1351.87854986277;1379.05991114751;314.771935700433;497.292290224317;-533.575868853380;345.335872655316;-3440.74107426101;1735.76314107882;-4242.63151155773].';
    % pole = 0.4; 
%     Qnum = [-582.2736;-116.7019;728.2703;-421.8695;343.5656;-919.1746;10.4039;-213.2097;-268.0179;-34.3649;99.3845;-237.2695;-575.8796;2246.0447;-1540.2489].';

    Lq = length(Qnum);

    %% Load the ID model from lib
    ID_Model;
    plant = plantX;

    %% Do Coprime Factorization  from lib
    plant = CoprimeFactorizationSS(plant,pole);

    %% Do Coprime Factorization via LFT from lib
    plant = LinearFractionalTransformation(plant);

    %% Check the constraints fitness
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
    % inputdata = myCosine(Test_time,Ts,Test_frequency,Test_amplitude);
    inputdata = mySine(Test_time,Ts,Test_frequency,Test_amplitude);

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

    % PlotTimeDomainResult('Fixed parameters - X',rk,yk,ek,uk,t);

%     max_ss = max(ek);
    max_mid_ss = max(ek(oneperiod));
%     rmse_ss = rms(ek);
    rmse_mid_ss = rms(ek(oneperiod));
    % fprintf("(stste space)\nL: %d\nMax error: %.6f um\nRMS(all): %.6f um\nRMS(onewave): %.6f um\n", Lq, max_ss, rmse_ss, rmse_mid_ss);

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
    
%     max_ss_all(idx) = max_ss;
    max_ss_one(idx) = max_mid_ss;
%     rmse_ss_all(idx) = rmse_ss;
    rmse_ss_one(idx) = rmse_mid_ss;
end

% figure();
% subplot(3,1,1);
% plot(max_ss_one,'-r','LineWidth',1.5); hold on; grid on; axis auto;
% ylabel('max (um)','FontSize',14);
% 
% % subplot(3,1,2);
% plot(rmse_ss_all,'-g','LineWidth',1.5); hold on; grid on; axis auto;
% ylabel('rmse (um)','FontSize',14);
% 
% % subplot(3,1,3);
% plot(rmse_mid_ss_all,':b','LineWidth',1.5); hold on; grid on; axis auto;
% ylabel('rmse period (um)','FontSize',14);

csvwrite('C:/Users/dalek669528/Desktop/CNC_proj/MATLAB/data/Generate/temp_error.csv', [max_ss_one';rmse_ss_one']');

% client = tcpclient('localhost', 50000);
write(client, uint8('k'))
exit;

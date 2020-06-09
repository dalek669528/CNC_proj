clear;clc;close all;
Ts = 0.001;
load('./ITRI_data/rawdata_1.mat')
load('./ITRI_data/rawdata_2.mat')
load('./ITRI_data/rawdata_3.mat')
load('./ITRI_data/rawdata_4.mat')

data1 = iddata(rawdata_1.X_Vout,rawdata_1.X_Vin,Ts);
data2 = iddata(rawdata_2.X_Vout,rawdata_2.X_Vin,Ts);
data3 = iddata(rawdata_3.X_Vout,rawdata_3.X_Vin,Ts);
data4 = iddata(rawdata_4.X_Vout,rawdata_4.X_Vin,Ts);
% iddata(y, u, Ts)
% creates an iddata object containing a time-domain output signal y and
% input signal u, respectively. Ts specifies the sample time of the
% experimental data

dataarray = [{data1},{data2},{data3},{data4}];
index = [1 2 3;1 2 4;1 3 4;2 3 4];
 
for delay = 0 %0:15
    for j = 1 %1:4
        idx = index(j,:);
        data_est = merge(dataarray{idx(1)},dataarray{idx(2)},dataarray{idx(3)});
%       m = merge(m1, m2, ...., mN)
%       merges the estimated models. The models m1,m2,...,mN must all be
%       the same structure, just diggering in parameter values and
%       covariance matrices.
%       m is the merged model, where the parameter vector is a
%       statistically weighted mean (using the covariance matrices to
%       determine the weights) of the parameter of mk

        for nx = 4 %3:1:5
%             sys = ssest(data_est,nx,...
%                         'Ts',Ts,...
%                         'InputDelay',delay,...
%                         'form','modal',...
%                         'DisturbanceModel','none'); 
            sys = tfest(data_est,nx,...
                        'Ts',Ts,...
                        'InputDelay',delay);
%           sys is a transfer function in z domain
                    
%           tfest(data, np, nz, ) : estimates a continuous-time transfer function, sys,
%           using time- or frequency-domain data, and contains np poles.
%           The number of zeros in the sys is max (np-1, 0)
                    
%           data : etimation data. For time domain estimation, data is an
%           iddata object containing the input and output signal values.
                    
                    
            if (max(abs(pole(sys)))<1 && max(abs(zero(sys)))<1)
                isminphase(j,nx,delay+1) = 1;
            else
                isminphase(j,nx,delay+1) = 0;
            end

            for i = 1:4
                fitness(i,nx,j,delay+1) = checkIDfitness(dataarray{i},sys);
            end
            fprintf("delay:%d idx:%d nx:%d\n",delay,j,nx);
        end
    end
end
opts1=bodeoptions('cstprefs');
opts1.XLim={[1e-02 3e03]};
bode(sys, opts1);grid on;
% margin(sys);grid on; % to see phase margin and gain margin
% kk = [fitness(:,3:5,1,b);isminphase(1,3:5,b);0 0 0;fitness(:,3:5,2,b);isminphase(2,3:5,b);0 0 0;fitness(:,3:5,3,b);isminphase(3,3:5,b);0 0 0;fitness(:,3:5,4,b);isminphase(4,3:5,b)];
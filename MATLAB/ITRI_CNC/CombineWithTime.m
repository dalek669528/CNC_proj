function [finaldata,T] = CombineWithTime(data,Ts)

t = 0:Ts:(length(data)-1)*Ts;
input = zeros(length(data),2);
input(:,1) = t;
input(:,2) = data;

finaldata = input;
T = (length(data)-1)*Ts;

end
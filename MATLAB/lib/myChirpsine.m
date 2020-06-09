function data = myChirpsine(T, Ts, f0, f1)

t = 0:Ts:T;
k = (f1-f0)/T;
phi = 2*pi*(f0+t*k/2).*t;
dphi = 2*pi*f0+2*k*pi*t+1;
A = 25./dphi;
% A = 1;
% A = 1./(1+t);
P =  A.*sin(phi);
input = zeros(length(t),2);
input(:,1) = t;
input(:,2) = P;

data = input;

end
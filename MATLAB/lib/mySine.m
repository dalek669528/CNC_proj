function data = mySine(T, Ts, w, amp)

t = 0:Ts:T;
P =  amp.*sin(w*t);
input = zeros(length(t),2);
input(:,1) = t;
input(:,2) = P;

data = input;

end
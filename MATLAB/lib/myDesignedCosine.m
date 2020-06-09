function data = myDesignedCosine(T, Ts, w, amp)

t = 0:Ts:T;
P =  0.5*amp.*(cos(w*t)+cos(2*w*t));
input = zeros(length(t),2);
input(:,1) = t;
input(:,2) = P-amp;

data = input;

end
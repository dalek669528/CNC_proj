%  LFT Expanded Discrete State Space 
%   Input: plant information, FIR coefficient
%   Output: Controller coefficient (State Space)
function [Ak,Bk,Ck,Dk] = LFTExpandedSS(J,Q)
    lq = length(Q);
    Q = reshape(Q,[1 lq]);
    qa1 = zeros(lq-1,1);
    qa2 = zeros(lq-1,1);
    qa1(2) = 1;
    QA = toeplitz(qa1,qa2);
    QB = [1;zeros(lq-2,1)]; 
    QC = Q(2:lq);
    QD = Q(1);
    
    F = 1/(1-QD*J.D(2,2));
    Dk = J.D(1,1) + J.D(1,2)*F*QD*J.D(2,1); 
    Ck = [J.C(1,:)+J.D(1,2)*F*QD*J.C(2,:),...
          J.D(1,2)*F*QC];
    Bk = [J.B(:,1)+J.B(:,2)*F*QD*J.D(2,1);...
          QB*J.D(2,1)+QB*J.D(2,2)*F*QD*J.D(2,1)];
    Ak = [J.A(:,:)+J.B(:,2)*F*QD*J.C(2,:),...
          J.B(:,2)*F*QC;...
          QB*J.C(2,:)+QB*J.D(2,2)*F*QD*J.C(2,:),...
          QA+QB*J.D(2,2)*F*QC];
end

         
                     
function [R]=transmatrixo(psiQ, thetaQ, phiQ)
%https://mathworld.wolfram.com/EulerAngles.html
%syms psiQ thetaQ phiQ
% syms x
% psiQ = 1.2;
% thetaQ = 1.3;
% phiQ = 1.1;
% x = 0.5;
method = 2;
switch method
    case 1
%% z-axis --> x-axis --> z-axis (x convention)
R1 = [cos(psiQ) -sin(psiQ) 0
    sin(psiQ) cos(psiQ) 0
    0 0 1];

R2 = [1 0 0
    0 cos(thetaQ) sin(thetaQ)
    0 -sin(thetaQ) cos(thetaQ)];

R3 = [cos(phiQ) sin(phiQ) 0
    -sin(phiQ) cos(phiQ) 0
    0 0 1];

    case 2

%% z-axis --> y-axis --> z-axis (x convention)
R1 = [cos(psiQ) sin(psiQ) 0
    -sin(psiQ) cos(psiQ) 0
    0 0 1];

R2 = [cos(thetaQ) 0 -sin(thetaQ)
    0 1 0
    sin(thetaQ) 0 cos(thetaQ)];

R3 = [cos(phiQ) sin(phiQ) 0
    -sin(phiQ) cos(phiQ) 0
    0 0 1];

    case 3
%% z-axis --> y-axis --> x-axis (x y z convention) pitch-roll-yaw
R1 = [cos(psiQ) -sin(psiQ) 0
    sin(psiQ) cos(psiQ) 0
    0 0 1];

R2 = [cos(thetaQ) 0 sin(thetaQ)
    0 1 0
    -sin(thetaQ) 0 cos(thetaQ)];

R3 = [1 0 0
    0 cos(phiQ) -sin(phiQ)
    0 sin(phiQ) cos(phiQ)];

end


R = R3*R2*R1;
end

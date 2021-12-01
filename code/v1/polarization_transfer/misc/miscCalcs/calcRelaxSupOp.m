function [G] = calcRelaxSupOp(T1_I,T2_I,T1_S,T2_S,Ix,Iy,Iz,Sx,Sy,Sz,nspins)
% Calculates relaxation superoperator for coupled spins in
% the uncerrolated random fields regime
%
% Higher order spin terms relax with the sum of hte individual rates:
%   R_IS = R_I+R_S
%%%%%%%%%%%%%%%%%%%%%%%

G = zeros((2^nspins)^2);
% calculating relaxation parameters
R1_I = 1/4/T1_I;
R1_S = 1/4/T1_S;
R2_I = 1/T2_I-1/2/T1_I;
R2_S = 1/T2_S-1/2/T1_S;
R3_I = 1/2/T1_I;
R3_S = 1/2/T1_S;

% constructing Zeeman basis
Sp = Sx + 1i*Sy;
Sn = Sx - 1i*Sy;
Ip = Ix + 1i*Iy;
In = Ix - 1i*Iy;


% S spins = protons
G = G + R2_S*kron(eye(2^nspins),Sz*Sz);
G = G + R1_S*kron(eye(2^nspins),Sp*Sn + Sn*Sp);
G = G + R3_S*kron(Sn.',Sp);
G = G - R3_S*kron(Sp.',Sn);
G = G + R3_S*kron(Sz.',eye(2^nspins));
G = G + R3_S*kron(eye(2^nspins),Sz);

% I spin = carbon
G = G + R2_I*kron(eye(2^nspins),Iz*Iz);
G = G + R1_I*kron(eye(2^nspins),Ip*In + In*Ip);
G = G + R3_I*kron(In.',Ip);
G = G - R3_I*kron(Ip.',In);
G = G + R3_I*kron(Iz.',eye(2^nspins));
G = G + R3_I*kron(eye(2^nspins),Iz);


end


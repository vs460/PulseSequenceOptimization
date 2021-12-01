function [L] = calcLiouvillian_formate(gamma_H,gamma_C,J,T1_I,T2_I,T1_S,T2_S,Ix,Iy,Iz,Sx,Sy,Sz,B1_L,B1_S,pulse_phase_L,pulse_phase_S,Omega_I,Omega_S)
% Caclulates Liouvilian for a simple IS spin system, e.g. 13C-formate

%% relaxation rates
%transverse
lambda_I = 1/T2_I;
lambda_S = 1/T2_S;
%longitudinal
rho_I = 1/T1_I;
rho_S = 1/T1_S;
%anti-phase
rho_Ia = lambda_I + rho_S;
rho_Sa = rho_I + lambda_S;
%multiple
lambda_mq = lambda_I + lambda_S;
rho_IS2sp = rho_I + rho_S;
%cross
sigma = 0;
mu_mq = 0;
delta_S = 0;
eta_S = 0;
%thermal correction
Theta_I = rho_I*0 + sigma*1;
Theta_S = rho_S*1 + sigma*0;
Theta_IS = delta_S*1;
%pulse
omega_Ix = 2*pi * gamma_C * B1_L * cos(pulse_phase_L);
omega_Iy = 2*pi * gamma_C * B1_L * sin(pulse_phase_L);
omega_Sx = 2*pi * gamma_H * B1_S * cos(pulse_phase_S);
omega_Sy = 2*pi * gamma_H * B1_S * sin(pulse_phase_S);


%% construct Liouvillian
L = [0              0         0        0      0        0         0        0         0       0          0         0            0       0           0      0;
     0          lambda_I   Omega_I -omega_Iy  0        0         0        0       pi*J      0          0         0            0       0           0      0;
     0          -Omega_I  lambda_I  omega_Ix  0        0         0      -pi*J       0       0          0         0            0       0           0      0;
     -2*Theta_I omega_Iy -omega_Ix   rho_I    0        0      sigma       0         0       0          0         0            0       0           0      0;
     0              0         0        0   lambda_S Omega_S -omega_Sy     0         0     eta_S       pi*J       0            0       0           0      0;
     0              0         0        0   -Omega_S lambda_S omega_Sx     0         0     -pi*J       eta_S      0            0       0           0      0;
     -2*Theta_S     0         0      sigma omega_Sy -omega_Sx  rho_S      0         0       0          0         0            0       0           0      delta_S;
     0              0       pi*J       0      0        0         0      rho_Ia   Omega_I    0          0     omega_Sy     -omega_Sx   0           0     -omega_Iy;
     0            -pi*J       0        0      0        0         0    -Omega_I    rho_Ia    0          0         0            0     omega_Sy -omega_Sx   omega_Ix;
     0              0         0        0     eta_S   pi*J        0        0         0     rho_Sa    Omega_S  omega_Iy         0    -omega_Ix      0     -omega_Sy;
     0              0         0        0     -pi*J   eta_S       0        0         0    -Omega_S    rho_Sa      0         omega_Iy   0      -omega_Ix   omega_Sx;
     0              0         0        0      0        0         0    -omega_Sy     0   -omega_Iy      0     lambda_mq     Omega_S  Omega_I    -mu_mq    0;
     0              0         0        0      0        0         0    omega_Sx      0        0     -omega_Iy -Omega_S     lambda_mq   mu_mq     Omega_I   0;
     0              0         0        0      0        0         0        0    -omega_Sy omega_Ix      0    -Omega_I        mu_mq   lambda_mq  Omega_S   0;
     0              0         0        0      0        0         0        0     omega_Sx     0     omega_Ix   -mu_mq       -Omega_I -Omega_S   lambda_mq 0;
     -2*Theta_IS    0         0        0      0        0      delta_S omega_Iy -omega_Ix omega_Sy -omega_Sx      0            0       0           0      rho_IS2sp];


end


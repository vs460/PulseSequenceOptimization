function[fitness] = polarizationTransfered_formate_relax_interleaved(population)
tic
[npop,nvar] = size(population);
gamma_C = 10.71e6;                 % 13C gyromagnetic ratio
gamma_H = 42.57e6;                 % 1H gyromagnetic ratio
omega_H = 600e6;                   % 1H resonance frequency for offresonance calculation
omega_C = omega_H/gamma_H*gamma_C; % 13C resonance frequency for offresonance calculation
J = 195;                           % coupling in Hz
B1_range = 0.3:0.1:1;                  % pulse amplitude range for compensation
fRange = -40:8:40;                   % offresonance range for compensation
T1_I = 15.8;%19
T2_I = 3.5;
T1_S = 2.2;%8
T2_S = 1.6;
% current desing works for I_nS_m spin systems
n_of_I_spins = 1;
n_of_S_spins = 1;
nspins = n_of_I_spins + n_of_S_spins;

% Define Pauli matrices
sigma_x = sparse([0 1/2; 1/2 0]);
sigma_y = sparse([0 -1i/2; 1i/2 0]);
sigma_z = sparse([1/2 0; 0 -1/2]);
unit = sparse([1 0; 0 1]);

% Calculate spin operators
Lx = zeros(2^nspins,2^nspins,nspins); Ly = zeros(2^nspins,2^nspins,nspins); Lz = zeros(2^nspins,2^nspins,nspins);
for l = 1:nspins
    Lx_current = 1;Ly_current = 1;Lz_current = 1;
    for s = 1:nspins
        if s == l
            Lx_current = kron(Lx_current,sigma_x);
            Ly_current = kron(Ly_current,sigma_y);
            Lz_current = kron(Lz_current,sigma_z);
        else
            Lx_current = kron(Lx_current,unit);
            Ly_current = kron(Ly_current,unit);
            Lz_current = kron(Lz_current,unit);
        end
    end
    Lx(:,:,l) = Lx_current; Ly(:,:,l) = Ly_current; Lz(:,:,l) = Lz_current;
end

% Use magnetic equivalence of the spins
Ix = sum(Lx(:,:,1:n_of_I_spins),3);
Iy = sum(Ly(:,:,1:n_of_I_spins),3);
Iz = sum(Lz(:,:,1:n_of_I_spins),3);

Sx = sum(Lx(:,:,n_of_I_spins+1:nspins),3);
Sy = sum(Ly(:,:,n_of_I_spins+1:nspins),3);
Sz = sum(Lz(:,:,n_of_I_spins+1:nspins),3);

% Detection state
coil = Ix + 1i*Iy;%Quadrature detection

% contaminating terms
other_terms = Sx+Sy+Sz+Ix*Sx+Ix*Sy+Ix*Sz+Iy*Sx+Iy*Sy+Iy*Sz+Iz*Sx+Iz*Sy+Iz*Sz;
%G_sup = (calcRelaxSupOp(T1_I,T2_I,T1_S,T2_S,Ix,Iy,Iz,Sx,Sy,Sz,nspins));

channelChange = 1e-3;
parfor i = 1:npop
    parameters = population(i,1:nvar);
    signal = zeros(length(B1_range),length(B1_range));
    signal_tmp = 0;
    for deltaf = fRange
        U = calcLiouvillian_formate(gamma_H,gamma_C,J,T1_I,T2_I,T1_S,T2_S,Ix,Iy,Iz,Sx,Sy,Sz,0,0,0,0,2*pi*deltaf,2*pi*deltaf*gamma_H/gamma_C);
        prop_free = expm(-U*channelChange);
        bb_C = 1;
        sum_weight = 0;
        for B1_C = B1_range
            bb_H = 1;
            for B1_H = B1_range
                % Initial state
                rho = zeros((2^nspins)^2,1);
                rho(7) = 1;
                % Pulse sequence
                for j = 1:6:nvar
                    % proton pulse
                    B1 = B1_H*parameters(j);
                    tau = parameters(j+2);
                    pulse_phase_H = parameters(j+1);
                    U = calcLiouvillian_formate(gamma_H,gamma_C,J,T1_I,T2_I,T1_S,T2_S,Ix,Iy,Iz,Sx,Sy,Sz,0,B1,0,pulse_phase_H,2*pi*deltaf,2*pi*deltaf*gamma_H/gamma_C);
                    rho = expm(-U*tau)*rho;
                    
                    %free precession
                    rho = prop_free*rho;
                    
                    % carbon pulse
                    j = j + 3;
                    B1 = B1_C*parameters(j);
                    tau = parameters(j+2);
                    pulse_phase_C = parameters(j+1);
                    U = calcLiouvillian_formate(gamma_H,gamma_C,J,T1_I,T2_I,T1_S,T2_S,Ix,Iy,Iz,Sx,Sy,Sz,B1,0,pulse_phase_C,0,2*pi*deltaf,2*pi*deltaf*gamma_H/gamma_C);
                    rho = expm(-U*tau)*rho;
                    
                    %free precession
                    rho = prop_free*rho;
                    
                end
                signal_tmp = (rho(2)+1i*rho(3));
                weight = 1/(1+(1-B1_C)^2+(1-B1_H)^2);
                signal(bb_C,bb_H) = signal(bb_C,bb_H) + weight*signal_tmp;
                sum_weight = sum_weight + weight;
                bb_H = bb_H + 1;
            end
            bb_C = bb_C +1;
        end
    end
    % normalization to 1 due to B1 and offres array
    signal = sum(sum(abs(signal)))/length(fRange)/sum_weight;
    fitness(i)	= 1e3*(-abs(signal));
end
toc
end







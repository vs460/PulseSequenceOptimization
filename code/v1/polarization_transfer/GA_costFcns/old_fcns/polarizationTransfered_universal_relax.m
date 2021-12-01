function[fitness] = polarizationTransfered_universal_relax(population)
tic
[npop,nvar] = size(population);
gamma_C = 10.71e6;                 % 13C gyromagnetic ratio
gamma_H = 42.57e6;                 % 1H gyromagnetic ratio
omega_H = 600e6;                   % 1H resonance frequency for offresonance calculation
omega_C = omega_H/gamma_H*gamma_C; % 13C resonance frequency for offresonance calculation
J = 195;                           % coupling in Hz
B1_range = 0.9:0.05:1.1;                  % pulse amplitude range for compensation
fRange = -5:1:5;                   % offresonance range for compensation
T1_I = 15.8;
T2_I = 3.5;
T1_S = 2.2;
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

% Define alternative basis
sigma_x = sparse([0 1; 0 0]);
sigma_y = sparse([0 0; 1 0]);
sigma_z = sparse([0 0; 0 1]);
unit = sparse([1 0; 0 0]);


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


parfor i = 1:npop
    parameters = population(i,1:nvar);
    signal = zeros(length(B1_range),length(B1_range));
    contamination = 0;
    bb_C = 1;
    for B1_C = B1_range
        bb_H = 1;
        for B1_H = B1_range
            signal_tmp = 0;
            for deltaf = fRange
                % Initial state
                rho = zeros((2^nspins)^2,1);
                rho(7) = 1;
                % Pulse sequence
                for j = 1:5:nvar
                    B1_L = B1_C*parameters(j+2);
                    pulse_phase_L = parameters(j+3);
                    B1_S = B1_H*parameters(j);
                    pulse_phase_S = parameters(j+1);
                    tau = parameters(j+4);
                    U = calcLiouvillian_formate(gamma_H,gamma_C,J,T1_I,T2_I,T1_S,T2_S,Ix,Iy,Iz,Sx,Sy,Sz,B1_L,B1_S,pulse_phase_L,pulse_phase_S,2*pi*deltaf,2*pi*deltaf*gamma_H/gamma_C);
                    rho = expm(-U*tau)*rho;
                end
                signal_tmp = signal_tmp + (rho(2)+1i*rho(3));
            end
            signal(bb_C,bb_H) = signal_tmp;
            bb_H = bb_H + 1;
        end
        bb_C = bb_C +1;
    end
    % normalization to 1 due to B1 and offres array
    signal = sum(sum(abs(signal)))/length(B1_range)^2/length(fRange);
    fitness(i)	= 1e3*(-abs(signal));
end
toc
end







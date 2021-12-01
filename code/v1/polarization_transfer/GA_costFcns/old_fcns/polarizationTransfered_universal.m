function[fitness] = polarizationTransfered_universal(population)
tic
[npop,nvar] = size(population);
gamma_C = 10.71e6;                 % 13C gyromagnetic ratio
gamma_H = 42.57e6;                 % 1H gyromagnetic ratio
omega_H = 600e6;                   % 1H resonance frequency for offresonance calculation
omega_C = omega_H/gamma_H*gamma_C; % 13C resonance frequency for offresonance calculation
J = 4.1;                           % coupling in Hz
B1_range = 1:1:1;                  % pulse amplitude range for compensation
fRange = -0:1:0;                   % offresonance range for compensation

% current desing works for I_nS_m spin systems
n_of_I_spins = 1;
n_of_S_spins = 3;
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


for i = 1:npop
    parameters = population(i,1:nvar);
    signal = 0;
    contamination = 0;
    for B1 = B1_range
        for deltaOmega = fRange
            % Initial state
            rho = Sz; %spins on the z-axis, starting polarization on proton
           
            % Pulse sequence
            for j = 1:5:nvar
                B1_L = B1*parameters(j+2);
                pulse_phase_L = parameters(j+3);
                B1_S = B1*parameters(j);
                pulse_phase_S = parameters(j+1);
                tau = parameters(j+4);
                
                H = (2*pi*gamma_C*B1_L*(cos(pulse_phase_L)*Ix + sin(pulse_phase_L)*Iy) +...              %pulse on carbon
                    2*pi*gamma_H*B1_S*(cos(pulse_phase_S)*Sx + sin(pulse_phase_S)*Sy) +...               %pulse on proton
                    2*pi*deltaOmega*(gamma_H/gamma_C*Sz + Iz) + ...                                      %off-resonance
                    2*pi*J*(Iz*Sz));                                                                     %coupling
                U = expm(-1i*H*tau);
                rho = U*rho*U';%U_conj;
            end
            
            % Detection
            signal = signal + trace(coil'*rho);
            contamination = contamination + trace(other_terms'*rho);
        end
    end
    % normalization to 1 due to B1 and offres array
    signal = signal/length(B1_range)/length(fRange);
    contamination = contamination/length(B1_range)/length(fRange);
    timeTot = sum(parameters(5:5:end)); % length of pulse sequence
    fitness(i)	= 1e3*(-real(signal));
end
toc
end






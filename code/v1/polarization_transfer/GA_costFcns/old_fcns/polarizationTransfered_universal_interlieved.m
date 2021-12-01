function[fitness] = polarizationTransfered_universal_interlieved(population)
tic
[npop,nvar] = size(population);
gamma_C = 10.71e6;                 % 13C gyromagnetic ratio
gamma_H = 42.57e6;                 % 1H gyromagnetic ratio
omega_H = 600e6;                   % 1H resonance frequency for offresonance calculation
omega_C = omega_H/gamma_H*gamma_C; % 13C resonance frequency for offresonance calculation
J = 4.5;                           % coupling in Hz
B1_range = 1:1:1;
fRange = -0:1:0;

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

% Use magnetic equivalence
Ix = sum(Lx(:,:,1:n_of_I_spins),3);
Iy = sum(Ly(:,:,1:n_of_I_spins),3);
Iz = sum(Lz(:,:,1:n_of_I_spins),3);

Sx = sum(Lx(:,:,n_of_I_spins+1:nspins),3);
Sy = sum(Ly(:,:,n_of_I_spins+1:nspins),3);
Sz = sum(Lz(:,:,n_of_I_spins+1:nspins),3);

% Detection state
coil = Sx + 1i*Sy;%Quadrature detection
% contaminating terms
other_terms = Ix+Iy+Iz+Ix*Sx+Ix*Sy+Ix*Sz+Iy*Sx+Iy*Sy+Iy*Sz+Iz*Sx+Iz*Sy+Iz*Sz;


parfor i = 1:npop
    parameters = population(i,1:nvar);
    signal = 0;
    contamination = 0;
    for B1_rel = B1_range
        for deltaFreq = fRange
            % Initial state
            rho = Iz;               %spins on the z-axis, starting polarization on carbon
            carbon = true;         %first pulse slot is on proton (see next if statement)
            
            % Pulse sequence
            for j = 1:3:nvar
                B1 = B1_rel*parameters(j);          
                pulse_phase = parameters(j+1);
                tau = parameters(j+2);
                if carbon
                    H = 2*pi*gamma_C*B1*(cos(pulse_phase)*Ix + sin(pulse_phase)*Iy) +...               %pulse on carbon
                        2*pi*deltaFreq*(gamma_H/gamma_C*Sz + Iz) + ...                                 %off-resonance
                        2*pi*J*(Iz*Sz);                                                                %coupling
                    U = expm(-1i*H*tau);                                                               % build propagator
                    carbon = false;                                                                    % change to proton channel
                else
                    H = 2*pi*gamma_H*B1*(cos(pulse_phase)*Sx + sin(pulse_phase)*Sy) +...               %pulse on proton
                        2*pi*deltaFreq*(gamma_H/gamma_C*Sz + Iz) + ...                                 %off-resonance
                        2*pi*J*(Iz*Sz);                                                                %coupling
                    U = expm(-1i*H*tau);                                                               % build propagator
                    carbon = true;                                                                     % change to proton channel
                end
                rho = U*rho*U';
                % free-precession during transmit channel change
                H_free = 2*pi*deltaFreq*(gamma_H/gamma_C*Sz + Iz) + ...                                %off-resonance
                         2*pi*J*(Iz*Sz);                                                               %coupling
                rho = expm(-1i*H_free*4e-6)*rho*expm(1i*H_free*4e-6);    
            end
            % Detection
            signal = signal + 1/(1+(1-B1_rel)^2+(deltaFreq/50)^2)*trace(coil'*rho); % quadratic weighting of the central regions
            contamination = contamination + trace(other_terms'*rho);
        end
    end
    % normalization to 1 due to B1 and offres array
    signal = signal/length(B1_range)/length(fRange);
    contamination = contamination/length(B1_range)/length(fRange);
    timeTot = sum(parameters(3:3:end)); % total time of the pulse sequence
    fitness(i)	= 1e3*(-real(signal)+abs(contamination));
end

toc
end






function[fitness] = polarizationTransfered(population)
tic
[npop,nvar] = size(population);
gamma_C = 10.71e6;%population(1,nvar-2); % 13C gyromagnetic ratio
gamma_H = 42.57e6;%population(1,nvar-1); % 13C gyromagnetic ratio
omega_H = 600e6;
omega_C = omega_H/gamma_H*gamma_C; 
J = 195;%population(1,nvar);
B1_range = 1:1:1;
fRange = -0:1:0;


% Define Pauli matrices
sigma_x = [0 1/2; 1/2 0];
sigma_y = [0 -1i/2; 1i/2 0];
sigma_z = [1/2 0; 0 -1/2];
unit = [1 0; 0 1];

% Calculate two-spin operators
Lx = kron(sigma_x,unit); Sx = kron(unit,sigma_x);
Ly = kron(sigma_y,unit); Sy = kron(unit,sigma_y);
Lz = kron(sigma_z,unit); Sz = kron(unit,sigma_z);

% Detection state
coil = Sx+1i*Sy;%Quadrature detection

parfor i = 1:npop
    parameters = population(i,1:nvar);
    signal = 0;
    
    for B1 = B1_range
        for deltaOmega = fRange
            
            % Initial state
            rho = Lz + 0*Sz; %spins on the z-axis
            % Pulse sequence
            for j = 1:5:nvar
                B1_L = B1*parameters(j);
                pulse_phase_L = parameters(j+1);
                B1_S = B1*parameters(j+2);
                pulse_phase_S = parameters(j+3);
                tau = parameters(j+4);
                
                H = 2*pi*gamma_H*B1_L*(cos(pulse_phase_L)*Lx + sin(pulse_phase_L)*Ly) +... %pulse on hydrogen
                    2*pi*gamma_C*B1_S*(cos(pulse_phase_S)*Sx + sin(pulse_phase_S)*Sy) +... %pulse on carbon
                    2*pi*deltaOmega*(Sz + gamma_H/gamma_C*Lz) + ...                                             %off-resonance
                    2*pi*J*(Lz*Sz);                                             %coupling
                U = expm(-1i*H*tau);
                %U_conj = expm(1i*H*tau);
                rho = U*rho*U';%U_conj;
            end
            % Detection
            signal = signal + trace(coil'*rho); % Detection = trace ... projection of the density on the detection state
        end
    end
    signal = signal/length(B1_range)/length(fRange);
    timeTot = sum(parameters(5:5:end));
    SAR = sum(parameters(1:5:end-4).*parameters(5:5:end)) + sum(parameters(3:5:end-2).*parameters(5:5:end)); %~ integrated power on both channels
    fitness(i)	= 1e3*(-abs(signal) + 100*timeTot);%abs(trace(((Lx+Ly+Lz+Sz)+2*(Lx*Sx+Ly*Sy+Lz*Sz+Lx*Sy+Lx*Sz+Ly*Sz+Lz*Sx+Ly*Sx))'*rho)));
    signalArray(i) = signal;
end
[M,I] =min(fitness);
signalArray(I)
toc
end




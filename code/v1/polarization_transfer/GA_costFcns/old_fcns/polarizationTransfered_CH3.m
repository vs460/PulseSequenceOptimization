function[fitness] = polarizationTransfered_CH3(population)
tic
[npop,nvar] = size(population);
gamma_C = 10.71e6;%population(1,nvar-2); % 13C gyromagnetic ratio
gamma_H = 42.57e6;%population(1,nvar-1); % 13C gyromagnetic ratio
omega_H = 600e6;
omega_C = omega_H/gamma_H*gamma_C; 
J = 4.1;%population(1,nvar);
B1_range = 0.5:0.2:2;
fRange = -50:10:50;


% Define Pauli matrices
sigma_x = sparse([0 1/2; 1/2 0]);
sigma_y = sparse([0 -1i/2; 1i/2 0]);
sigma_z = sparse([1/2 0; 0 -1/2]);
unit = sparse([1 0; 0 1]);

% Calculate four-spin operators
%carbon
Lx = kron(kron(kron(sigma_x,unit),unit),unit); 
Ly = kron(kron(kron(sigma_y,unit),unit),unit);  
Lz = kron(kron(kron(sigma_z,unit),unit),unit);
%proton 1
Sx_1 = kron(kron(kron(unit,sigma_x),unit),unit); 
Sy_1 = kron(kron(kron(unit,sigma_y),unit),unit);
Sz_1 = kron(kron(kron(unit,sigma_z),unit),unit);
%proton 2
Sx_2 = kron(kron(unit,kron(unit,sigma_x)),unit); 
Sy_2 = kron(kron(unit,kron(unit,sigma_y)),unit); 
Sz_2 = kron(kron(unit,kron(unit,sigma_z)),unit); 
%proton3 
Sx_3 = kron(unit,kron(unit,kron(unit,sigma_x)));
Sy_3 = kron(unit,kron(unit,kron(unit,sigma_y)));
Sz_3 = kron(unit,kron(unit,kron(unit,sigma_z)));

Sx = (Sx_1+Sx_2+Sx_3);
Sy = (Sy_1+Sy_2+Sy_3);
Sz = (Sz_1+Sz_2+Sz_3);

% Detection state
coil = Sx + 1i*Sy;%Quadrature detection

%contaminating states
other_terms = ((Lx+Ly+Lz)+...
   2*(Lx*Sx+Ly*Sy+Lz*Sz...
   +Lx*Sy+Lx*Sz+Ly*Sz+Ly*Sx+Lz*Sx+Lz*Sy));


parfor i = 1:npop
    parameters = population(i,1:nvar);
    signal = 0;
    contamination = 0;
    for B1 = B1_range
        for deltaOmega = fRange
            % Initial state
            rho = Lz; %spins on the z-axis, starting polarization on carbon
           
            % Pulse sequence
            for j = 1:5:nvar
                B1_L = B1*parameters(j+2);
                pulse_phase_L = parameters(j+3);
                B1_S = B1*parameters(j);
                pulse_phase_S = parameters(j+1);
                tau = parameters(j+4);
                
                H = (2*pi*gamma_C*B1_L*(cos(pulse_phase_L)*Lx + sin(pulse_phase_L)*Ly) +...              %pulse on carbon
                    2*pi*gamma_H*B1_S*(cos(pulse_phase_S)*Sx + sin(pulse_phase_S)*Sy) +...               %pulse on proton
                    2*pi*deltaOmega*(gamma_H/gamma_C*Sz + Lz) + ...                                      %off-resonance
                    2*pi*J*(Lz*Sz));                                                                     %coupling
                U = expm(-1i*H*tau);
                rho = U*rho*U';%U_conj;
            end
            
            % Detection
            signal = signal + coil'*rho; % Detection = trace ... projection of the density on the detection state
            contamination = contamination + other_terms'*rho;
        end
    end
    signal = signal/length(B1_range)/length(fRange);
    contamination = contamination/length(B1_range)/length(fRange);
    timeTot = sum(parameters(5:5:end));
    %SAR = sum(parameters(1:5:end-4).*parameters(5:5:end)) + sum(parameters(3:5:end-2).*parameters(5:5:end)); %~ integrated power on both channels
    fitness(i)	= 1e3*(-abs(trace(signal))+timeTot);
end
toc
end





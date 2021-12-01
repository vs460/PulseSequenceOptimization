function[fitness] = polarizationTransfered_universal_relax_diffBasis(population)
tic

%% load pre-generated Hamiltonian superoperators for the faseter computation
% only those components are needed that are present in the
% Hamiltonian...taking advantage of linearity
genH_sup = false;
if ~genH_sup
    path = 'C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\NaAcetate';
    load(strcat(path,'\H_sup_Ix.mat'));
    load(strcat(path,'\H_sup_Iy.mat'));
    load(strcat(path,'\H_sup_Iz.mat'));
    load(strcat(path,'\H_sup_Sx.mat'));
    load(strcat(path,'\H_sup_Sy.mat'));
    load(strcat(path,'\H_sup_Sz.mat'));
    load(strcat(path,'\H_sup_IzSz.mat'));
    H_sup_Ix = (double(H_sup_Ix));
    H_sup_Iy = (double(H_sup_Iy));
    H_sup_Iz = (double(H_sup_Iz));
    H_sup_Sx = (double(H_sup_Sx));
    H_sup_Sy = (double(H_sup_Sy));
    H_sup_Sz = (double(H_sup_Sz));
    H_sup_IzSz = (double(H_sup_IzSz));
end

[npop,nvar] = size(population);
fitness = zeros(npop,1);
gamma_C = 10.71e6;                 % 13C gyromagnetic ratio
gamma_H = 42.57e6;                 % 1H gyromagnetic ratio
omega_H = 600e6;                   % 1H resonance frequency for offresonance calculation
omega_C = omega_H/gamma_H*gamma_C; % 13C resonance frequency for offresonance calculation
J = 6;                           % coupling in Hz
B1_range = 0.5:0.25:1.5;                  % pulse amplitude range for compensation
fRange = -40:20:40;                   % offresonance range for compensation
channelChange = 400e-6;

% relaxation parameters
T1_I = 15.8;                       % 13C longitudinal
T2_I = 3.5;                        % 13C transverse
T1_S = 2.2;                        % 1H longitudinal
T2_S = 1.6;                        % 1H transverse
% relaxation rates
R1_I = 1/T1_I; R2_I = 1/T2_I; R1_S = 1/T1_S; R2_S = 1/T2_S;

% **********************************************************
% ******* current desing works for I1S3 spin systems *******
% **********************************************************
n_of_I_spins = 1;
n_of_S_spins = 3;
nspins = n_of_I_spins + n_of_S_spins;

% relaxation rates stored in a matrix for the product operator rate
% calculations
relaxation_rates = repmat([0 R2_I R2_I R1_I],[n_of_I_spins,1]);
relaxation_rates = [relaxation_rates;repmat([0 R2_S R2_S R1_S],[n_of_S_spins,1])];

% Define Pauli matrices
sigma_x =  ([0 1/2; 1/2 0]);
sigma_y =  ([0 -1i/2; 1i/2 0]);
sigma_z =  ([1/2 0; 0 -1/2]);
unit =  ([1 0; 0 1]);

% Calculate spin operators
L = zeros(2^nspins,2^nspins,nspins,4);
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
    L(:,:,l,1) = eye(2^nspins); L(:,:,l,2) = Lx_current; L(:,:,l,3) = Ly_current; L(:,:,l,4) = Lz_current;
end

% Use magnetic equivalence of the spins
Ix = sum(L(:,:,1:n_of_I_spins,2),3);
Iy = sum(L(:,:,1:n_of_I_spins,3),3);
Iz = sum(L(:,:,1:n_of_I_spins,4),3);

Sx = sum(L(:,:,n_of_I_spins+1:nspins,2),3);
Sy = sum(L(:,:,n_of_I_spins+1:nspins,3),3);
Sz = sum(L(:,:,n_of_I_spins+1:nspins,4),3);

% Detection state
coil = Ix + 1i*Iy;%Quadrature detection

% construct cartesian product operator basis
% basis elements are represented as a 1x4 row vector:
% entry 1 = E, 2 = Lx, 3 = Ly, 4 = Lz
% and the basis operators are the product of the permutations (with repetitions) of these
for bb = 1:nspins
    permutation_mtx_tmp(:,bb) = reshape(repmat(mod((1:4^(nspins-bb+1))-1,4)'+1,[1,4^(bb-1)])',[1,4^nspins])';
end

% let's order the basis operators so that E and Ixyz and S123xyz are the
% first 13 entries
permutation_mtx = zeros(size(permutation_mtx_tmp));
permutation_mtx_tmp(1:4,:) = [];            %removing Ix Iy Iz
permutation_mtx_tmp(5-4,:) = [];            % S1x ....substracted number = #of already deleted rows
permutation_mtx_tmp(5+4-5,:) = [];          % S1y
permutation_mtx_tmp(5+2*4-6,:) = [];        % S1z
permutation_mtx_tmp(17-7,:) = [];           % S2x
permutation_mtx_tmp(17+4^2-8,:) = [];       % S2y
permutation_mtx_tmp(17+2*4^2-9,:) = [];     % S2z
permutation_mtx_tmp(4^3+1-10,:) = [];       % S1x
permutation_mtx_tmp(4^3+1+4^3-11,:) = [];   % S1y
permutation_mtx_tmp(4^3+1+2*4^3-12,:) = []; % S1z

% write the E, Ixyz and S123yz operators to the first 13 rows
permutation_mtx(1,:) = [1;1;1;1];
permutation_mtx(2:4,:) = [2,1,1,1;3,1,1,1;4,1,1,1];
permutation_mtx(5:7,:) = circshift([2,1,1,1;3,1,1,1;4,1,1,1],[0,1]);
permutation_mtx(8:10,:) = circshift([2,1,1,1;3,1,1,1;4,1,1,1],[0,2]);
permutation_mtx(11:13,:) = circshift([2,1,1,1;3,1,1,1;4,1,1,1],[0,3]);
permutation_mtx(14:end,:) = permutation_mtx_tmp;

% constructing cell with the basis operators and the relaxation supermatrix
B = cell((2^nspins)^2,1);
R = zeros((2^nspins)^2,1);
B_block = [];
for bb = 1:4^nspins
    basis_prod_tmp = eye(2^nspins);
    for ns = 1:nspins
        basis_prod_tmp = basis_prod_tmp*L(:,:,ns,permutation_mtx(bb,ns)); % basis operator is the product of 4 terms according to the correspopnding permutatin
        R(bb) = R(bb) + relaxation_rates(ns,permutation_mtx(bb,ns));      % relaxation rate of product operators = sum of correspoding single rates
    end
    if bb >13    % normalization according to spindynamica exported basis set
        B(bb) = {basis_prod_tmp*2^(sum(permutation_mtx(bb,:)>1)-2)};
    elseif bb > 1
        B(bb) = {basis_prod_tmp/2};
    else
        B(bb) = {basis_prod_tmp/4};
    end
end
%B = cell((2^nspins)^2,1);
%B(1) = {eye(4)/2}; B(2) = {Ix}; B(3) = {Iy}; B(4) = {Iz}; B(5) = {Sx}; B(6) = {Sy}; B(7) = {Sz}; B(8) = {2*Ix*Sx};
%B(9) = {2*Ix*Sy}; B(10) = {2*Ix*Sz}; B(11) = {2*Iy*Sx}; B(12) = {2*Iy*Sy}; B(13) = {2*Iy*Sz}; B(14) = {2*Iz*Sx}; B(15) = {2*Iz*Sy}; B(16) = {2*Iz*Sz};
%R = [0 R2_I R2_I R1_I R2_S R2_S R1_S R2_I+R2_S R2_I+R2_S R2_I+R1_S R2_I+R2_S R2_I+R2_S R2_I+R1_S R1_I+R2_S R1_I+R2_S R1_I+R1_S];
% high temperature thermal correction
R = diag(R);
R(4,1) = -2*R1_I;
R(7,1) = -2*R1_S;
R(10,1) = -2*R1_S;
R(13,1) = -2*R1_S;

%% generating Hamiltonian supermatrix structure
if genH_sup
    H_Ix = Ix;
    H_Iy = Iy;
    H_Iz = Iz;
    H_Sx = Sx;
    H_Sy = Sy;
    H_Sz = Sz;
    H_IzSz = Iz*Sz;
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Ix(r,s) = trace(cell2mat(B(r))'*(H_Ix*cell2mat(B(s))-cell2mat(B(s))*H_Ix));
        end
    end
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Iy(r,s) = trace(cell2mat(B(r))'*(H_Iy*cell2mat(B(s))-cell2mat(B(s))*H_Iy));
        end
    end
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Iz(r,s) = trace(cell2mat(B(r))'*(H_Iz*cell2mat(B(s))-cell2mat(B(s))*H_Iz));
        end
    end
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Sx(r,s) = trace(cell2mat(B(r))'*(H_Sx*cell2mat(B(s))-cell2mat(B(s))*H_Sx));
        end
    end
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Sy(r,s) = trace(cell2mat(B(r))'*(H_Sy*cell2mat(B(s))-cell2mat(B(s))*H_Sy));
        end
    end
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Sz(r,s) = trace(cell2mat(B(r))'*(H_Sz*cell2mat(B(s))-cell2mat(B(s))*H_Sz));
        end
    end
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_IzSz(r,s) = trace(cell2mat(B(r))'*(H_IzSz*cell2mat(B(s))-cell2mat(B(s))*H_IzSz));
        end
    end
   save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\NaAcetate\H_sup_Ix','H_sup_Ix')
   save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\NaAcetate\H_sup_Iy','H_sup_Iy')
   save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\NaAcetate\H_sup_Iz','H_sup_Iz')
   save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\NaAcetate\H_sup_Sx','H_sup_Sx')
   save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\NaAcetate\H_sup_Sy','H_sup_Sy')
   save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\NaAcetate\H_sup_Sz','H_sup_Sz')
   save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\NaAcetate\H_sup_IzSz','H_sup_IzSz')
end
H_offres_sup = (H_sup_Iz + gamma_H/gamma_C*H_sup_Sz);
%% fitness evaluation
parfor i = 1:npop
    parameters = population(i,1:nvar);
    signal = zeros(length(B1_range),length(B1_range));
    contamination = 0;
    cos_pulse_phase = cos(parameters((1:3:nvar)+1));
    sin_pulse_phase = sin(parameters((1:3:nvar)+1));
    signal_tmp = 0;
    for deltaf = fRange
        H_sup_free = 2*pi*deltaf*H_offres_sup + ...                                                   %frequency offset
            2*pi*J*(H_sup_IzSz);                                                                      %coupling
        prop_free = (expm(((-1i*H_sup_free - R)*channelChange)));
        bb_C = 1;
        sum_weight = 0;
        for B1_C = B1_range
            bb_H = 1;
            for B1_H = B1_range
                % Initial state = spins along z-axis, polariation on protons
                rho = zeros((2^nspins)^2,1);
                rho(7) = 1;
                rho(10) = 1;
                rho(13) = 1;
                carbon = false;
                k = 1;
                % Pulse sequence
                for j = 1:6:nvar
                    B1 = B1_H*parameters(j);
                    tau = parameters(j+2);
                    
                    % proton
                    H_sup = (2*pi*gamma_H*B1*(cos_pulse_phase(k)*H_sup_Sx + sin_pulse_phase(k)*H_sup_Sy) +...       %pulse on carbon
                        2*pi*deltaf*H_offres_sup + ...                                                              %frequency offset
                        2*pi*J*(H_sup_IzSz));                                                                       %coupling

                    rho = expm(((-1i*H_sup - R)*tau))*rho;                                                            %build propagator
                    %exp_sup = expm(((-1i*H_sup - R)*tau));
                    %sum([((exp_sup(:,2)>0) + (prop_free(:,2)>0));(exp_sup(:,3)>0 + prop_free(:,3)>0)]>0)
                    % free-precession during transmit channel change
                    rho = prop_free*rho;
                    
                    % carbon
                    j = j + 3;
                    k = k + 1;
                    B1 = B1_C*parameters(j);
                    tau = parameters(j+2);
                    H_sup = (2*pi*gamma_C*B1*(cos_pulse_phase(k)*H_sup_Ix + sin_pulse_phase(k)*H_sup_Iy) +...     %pulse on carbon
                        2*pi*deltaf*H_offres_sup + ...                                                            %frequency offset
                        2*pi*J*(H_sup_IzSz));                                                                     %coupling

                    rho = expm(((-1i*H_sup - R)*tau))*rho;                                                          %build propagator
                    
                    % free-precession during transmit channel change
                    rho = prop_free*rho;
                    k = k + 1;
                    %exp_sup = expm(((-1i*H_sup - R)*tau));
                    %sum([((exp_sup(:,2)>0) + (prop_free(:,2)>0));(exp_sup(:,3)>0 + prop_free(:,3)>0)]>0)
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
    fitness(i)	= 1e3*(-real(signal)+contamination);
end
toc
end








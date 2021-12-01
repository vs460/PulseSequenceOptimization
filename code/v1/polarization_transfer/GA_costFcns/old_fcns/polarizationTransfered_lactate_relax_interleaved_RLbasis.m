function[fitness] = polarizationTransfered_lactate_relax_interleaved_RLbasis(population)
tic

%% load pre-generated Hamiltonian superoperators for the faseter computation
% only those components are needed that are present in the
% Hamiltonian...taking advantage of linearity
genH_sup = false;
if ~genH_sup
    path = 'C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization';
    load(fullfile(path,'H_sup_Im.mat'));
    load(fullfile(path,'H_sup_Ip.mat'));
    load(fullfile(path,'H_sup_Iz.mat'));
    load(fullfile(path,'H_sup_Sm_CH.mat'));
    load(fullfile(path,'H_sup_Sp_CH.mat'));
    load(fullfile(path,'H_sup_Sz_CH.mat'));
    load(fullfile(path,'H_sup_Sm_CH3.mat'));
    load(fullfile(path,'H_sup_Sp_CH3.mat'));
    load(fullfile(path,'H_sup_Sz_CH3.mat'));
    load(fullfile(path,'H_sup_IzSz_CH3.mat'));
    load(fullfile(path,'H_sup_IzSz_CH.mat'));
    load(fullfile(path,'H_sup_Sz_CHSz_CH3.mat'));
    
    H_sup_Im = double(H_sup_Im);
    H_sup_Ip = double(H_sup_Ip);
    H_sup_Iz = double(H_sup_Iz);
    H_sup_Sm_CH = double(H_sup_Sm_CH);
    H_sup_Sp_CH = double(H_sup_Sp_CH);
    H_sup_Sz_CH = double(H_sup_Sz_CH);
    H_sup_Sm_CH3 = double(H_sup_Sm_CH3);
    H_sup_Sp_CH3 = double(H_sup_Sp_CH3);
    H_sup_Sz_CH3 = double(H_sup_Sz_CH3);
    H_sup_IzSz_CH = double(H_sup_IzSz_CH);
    H_sup_IzSz_CH3 = double(H_sup_IzSz_CH3);
    H_sup_Sz_CHSz_CH3 = double(H_sup_Sz_CHSz_CH3);
end

[npop,nvar] = size(population);
fitness = zeros(npop,1);
gamma_C = 10.71e6;                 % 13C gyromagnetic ratio
gamma_H = 42.57e6;                 % 1H gyromagnetic ratio
omega_H = 600e6;                   % 1H resonance frequency for offresonance calculation
omega_C = omega_H/gamma_H*gamma_C; % 13C resonance frequency for offresonance calculation
J_CH3 = 4.1;                           % coupling in Hz
J_CH = 3.1;                           % coupling in Hz
J_HH3 = 7;                           % coupling in Hz
B1_range = 0.5:0.25:1.5;                  % pulse amplitude range for compensation
fRange = -40:40:40;                   % offresonance range for compensation
channelChange = 1e-3;
methinOffset = 2.6*600;

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
n_of_S_spins = 4;
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
    Lm_current = 1;Lp_current = 1;Lz_current = 1;
    for s = 1:nspins
        if s == l
            Lm_current = kron(Lm_current,sigma_x);
            Lp_current = kron(Lp_current,sigma_y);
            Lz_current = kron(Lz_current,sigma_z);
        else
            Lm_current = kron(Lm_current,unit);
            Lp_current = kron(Lp_current,unit);
            Lz_current = kron(Lz_current,unit);
        end
    end
    L(:,:,l,1) = eye(2^nspins); L(:,:,l,2) = Lm_current; L(:,:,l,3) = Lp_current; L(:,:,l,4) = Lz_current;
end

% Use magnetic equivalence of the spins
Im = sum(L(:,:,1:n_of_I_spins,2),3);
Ip = sum(L(:,:,1:n_of_I_spins,3),3);
Iz = sum(L(:,:,1:n_of_I_spins,4),3);

Sm_CH3 = sum(L(:,:,n_of_I_spins+1:nspins-1,2),3);
Sp_CH3 = sum(L(:,:,n_of_I_spins+1:nspins-1,3),3);
Sz_CH3 = sum(L(:,:,n_of_I_spins+1:nspins-1,4),3);

Sm_CH = L(:,:,nspins,2);
Sp_CH = L(:,:,nspins,3);
Sz_CH = L(:,:,nspins,4);

% construct cartesian product operator basis
% basis elements are represented as a 1x5 row vector:
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
permutation_mtx_tmp(4^3+1-10,:) = [];       % S3x
permutation_mtx_tmp(4^3+1+4^3-11,:) = [];   % S3y
permutation_mtx_tmp(4^3+1+2*4^3-12,:) = []; % S3z
permutation_mtx_tmp(4^4+1-13,:) = [];       % S3x_methin
permutation_mtx_tmp(4^4+1+4^4-14,:) = [];   % S3y_methin
permutation_mtx_tmp(4^4+1+2*4^4-15,:) = []; % S3z_methin


% write the E, Ixyz and S123yz operators to the first 13 rows
permutation_mtx(1,:) = [1,1,1,1,1];
permutation_mtx(2:4,:) = [2,1,1,1,1;3,1,1,1,1;4,1,1,1,1];
permutation_mtx(5:7,:) = circshift([2,1,1,1,1;3,1,1,1,1;4,1,1,1,1],[0,1]);
permutation_mtx(8:10,:) = circshift([2,1,1,1,1;3,1,1,1,1;4,1,1,1,1],[0,2]);
permutation_mtx(11:13,:) = circshift([2,1,1,1,1;3,1,1,1,1;4,1,1,1,1],[0,3]);
permutation_mtx(14:16,:) = circshift([2,1,1,1,1;3,1,1,1,1;4,1,1,1,1],[0,4]);
permutation_mtx(17:end,:) = permutation_mtx_tmp;

% constructing cell with the basis operators and the relaxation supermatrix
B = cell((2^nspins)^2,1);
R = zeros((2^nspins)^2,1);
B_block = [];
for bb = 1:4^nspins
    basis_prod_tmp = eye(2^nspins);
    for ns = 1:nspins
        basis_prod_tmp = basis_prod_tmp*L(:,:,ns,permutation_mtx(bb,ns)); % basis operator is the product of 5 terms according to the correspopnding permutation
        R(bb) = R(bb) + relaxation_rates(ns,permutation_mtx(bb,ns));      % relaxation rate of product operators = sum of correspoding single rates
    end
    if bb >16    % normalization according to spindynamica exported basis set
        B(bb) = {basis_prod_tmp*2^(sum(permutation_mtx(bb,:)>1)-2)/sqrt(2)};
    elseif bb > 1
        B(bb) = {basis_prod_tmp/2/sqrt(2)};
    else
        B(bb) = {basis_prod_tmp/4/sqrt(2)};
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
R(16,1) = -2*R1_S;


%% generating Hamiltonian supermatrix structure
if genH_sup
    h = waitbar(0,'Calculating superoperators...');
    H_Im = Im;
    H_Ip = Ip;
    H_Iz = Iz;
    H_Sm_CH3 = Sm_CH3;
    H_Sp_CH3 = Sp_CH3;
    H_Sz_CH3 = Sz_CH3;
    H_Sm_CH = Sm_CH;
    H_Sp_CH = Sp_CH;
    H_Sz_CH = Sz_CH;
    H_IzSz_CH3 = Iz*Sz_CH3;
    H_IzSz_CH = Iz*Sz_CH;
    H_Sz_CHSz_CH3 = Sz_CH*Sz_CH3;
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Im(r,s) = trace(cell2mat(B(r))'*(H_Im*cell2mat(B(s))-cell2mat(B(s))*H_Im));
        end
    end
    save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\H_sup_Im','H_sup_Im')
    waitbar(1/12,h,['Calculating superoperators: ' sprintf('%i%% along',round(1/12*100))]);
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Ip(r,s) = trace(cell2mat(B(r))'*(H_Ip*cell2mat(B(s))-cell2mat(B(s))*H_Ip));
        end
    end
    save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\H_sup_Ip','H_sup_Ip')
    waitbar(2/12,h,['Calculating superoperators: ' sprintf('%i%% along',round(2/12*100))]);
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Iz(r,s) = trace(cell2mat(B(r))'*(H_Iz*cell2mat(B(s))-cell2mat(B(s))*H_Iz));
        end
    end
    save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\H_sup_Iz','H_sup_Iz')
    waitbar(3/12,h,['Calculating superoperators: ' sprintf('%i%% along',round(3/12*100))]);
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Sm_CH3(r,s) = trace(cell2mat(B(r))'*(H_Sm_CH3*cell2mat(B(s))-cell2mat(B(s))*H_Sm_CH3));
        end
    end
    save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\H_sup_Sm_CH3','H_sup_Sm_CH3')
    waitbar(4/12,h,['Calculating superoperators: ' sprintf('%i%% along',round(4/12*100))]);
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Sp_CH3(r,s) = trace(cell2mat(B(r))'*(H_Sp_CH3*cell2mat(B(s))-cell2mat(B(s))*H_Sp_CH3));
        end
    end
    save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\H_sup_Sp_CH3','H_sup_Sp_CH3')
    waitbar(5/12,h,['Calculating superoperators: ' sprintf('%i%% along',round(5/12*100))]);
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Sz_CH3(r,s) = trace(cell2mat(B(r))'*(H_Sz_CH3*cell2mat(B(s))-cell2mat(B(s))*H_Sz_CH3));
        end
    end
    save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\H_sup_Sz_CH3','H_sup_Sz_CH3')
    waitbar(6/12,h,['Calculating superoperators: ' sprintf('%i%% along',round(6/12*100))]);
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_IzSz_CH3(r,s) = trace(cell2mat(B(r))'*(H_IzSz_CH3*cell2mat(B(s))-cell2mat(B(s))*H_IzSz_CH3));
        end
    end
    save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\H_sup_IzSz_CH3','H_sup_IzSz_CH3')
    waitbar(7/12,h,['Calculating superoperators: ' sprintf('%i%% along',round(7/12*100))]);
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_IzSz_CH(r,s) = trace(cell2mat(B(r))'*(H_IzSz_CH*cell2mat(B(s))-cell2mat(B(s))*H_IzSz_CH));
        end
    end
    save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\H_sup_IzSz_CH','H_sup_IzSz_CH')
    waitbar(8/12,h,['Calculating superoperators: ' sprintf('%i%% along',round(8/12*100))]);
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Sz_CHSz_CH3(r,s) = trace(cell2mat(B(r))'*(H_Sz_CHSz_CH3*cell2mat(B(s))-cell2mat(B(s))*H_Sz_CHSz_CH3));
        end
    end
    save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\H_sup_Sz_CHSz_CH3','H_sup_Sz_CHSz_CH3')
    waitbar(9/12,h,['Calculating superoperators: ' sprintf('%i%% along',round(9/12*100))]);
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Sm_CH(r,s) = trace(cell2mat(B(r))'*(H_Sm_CH*cell2mat(B(s))-cell2mat(B(s))*H_Sm_CH));
        end
    end
    save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\H_sup_Sm_CH','H_sup_Sm_CH')
    waitbar(10/12,h,['Calculating superoperators: ' sprintf('%i%% along',round(10/12*100))]);
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Sp_CH(r,s) = trace(cell2mat(B(r))'*(H_Sp_CH*cell2mat(B(s))-cell2mat(B(s))*H_Sp_CH));
        end
    end
    save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\H_sup_Sp_CH','H_sup_Sp_CH')
    waitbar(11/12,h,['Calculating superoperators: ' sprintf('%i%% along',round(11/12*100))]);
    for r = 1:(2^nspins)^2
        for s = 1:(2^nspins)^2
            H_sup_Sz_CH(r,s) = trace(cell2mat(B(r))'*(H_Sz_CH*cell2mat(B(s))-cell2mat(B(s))*H_Sz_CH));
        end
    end
    save('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization\H_sup_Sz_CH','H_sup_Sz_CH')
    waitbar(12/12,h,['Calculating superoperators: ' sprintf('%i%% along',round(12/12*100))]);
    close(h)
end
H_sup_Sm = (H_sup_Sm_CH+H_sup_Sm_CH3);
H_sup_Sp = (H_sup_Sp_CH+H_sup_Sp_CH3);
H_sup_Sz = (H_sup_Sz_CH+H_sup_Sz_CH3);
H_sup_Jcoupling = (J_CH3*H_sup_IzSz_CH3+J_CH*H_sup_IzSz_CH+J_HH3*H_sup_Sz_CHSz_CH3);
H_sup_offres = (H_sup_Iz + gamma_H/gamma_C*(H_sup_Sz_CH+H_sup_Sz_CH3));

%% fitness evaluation
for i = 1:npop
    parameters = population(i,1:nvar);
    signal = zeros(length(B1_range),length(B1_range));
    contamination = 0;
    cos_pulse_phase = cos(parameters((1:3:nvar)+1));
    sin_pulse_phase = sin(parameters((1:3:nvar)+1));
    signal_tmp = 0;
    for deltaf = fRange
        H_sup_free = 2*pi*deltaf*H_sup_offres + 2*pi*methinOffset*H_sup_Sz_CH + ...                %frequency offset
            2*pi*H_sup_Jcoupling;                                                                %coupling
        prop_free = (expm(gpuArray((-1i*H_sup_free - R)*channelChange)));
        bb_C = 1;
        sum_weight = 0;
        for B1_C = B1_range
            bb_H = 1;
            for B1_H = B1_range
                % Initial state = spins along z-axis, polariation on protons
                rho = zeros((2^nspins)^2,1,'gpuArray');
                rho(7) = 1;
                rho(10) = 1;
                rho(13) = 1;
                rho(16) = 1;
                carbon = false;
                k = 1;
                % Pulse sequence
                for j = 1:6:nvar
                    B1 = B1_H*parameters(j);
                    tau = parameters(j+2);
                    
                    % proton
                    H_sup = 2*pi*gamma_H*B1*(cos_pulse_phase(k)*(H_sup_Sp+H_sup_Sm)/sqrt(2) + sin_pulse_phase(k)*(H_sup_Sp-H_sup_Sm)/1i/sqrt(2)) + ...  %pulse on proton
                        2*pi*deltaf*H_sup_offres + 2*pi*methinOffset*H_sup_Sz_CH + ...                         %frequency offset
                        2*pi*H_sup_Jcoupling;                                                                %coupling
                    if (j == 1 && deltaf == min(fRange))
                        M_H = perm_mat_to_make_block_diag(H_sup);
                        R1 = R(1:512,1:512);
                        R2 = R(513:end,513:end);
                    end
                    block_H = gpuArray(M_H*(-1i*H_sup - R)*M_H');
                    block1_H = block_H(1:512,1:512);
                    block2_H = block_H(513:end,513:end);
                    block1_H_exp = (expm((block1_H*tau)));
                    block2_H_exp = (expm((block2_H*tau)));
                    prop_H = blkdiag(block1_H_exp,block2_H_exp);
                    rho = (M_H'*prop_H*M_H)*rho;                                                       %build propagator
                    
                    % free-precession during transmit channel change
                    rho = prop_free*rho;                                        
                    
                    % carbon
                    
                    k = k + 1;
                    B1 = B1_C*parameters(j+3);
                    tau = parameters(j+2+3);
                    H_sup = 2*pi*gamma_C*B1*(cos_pulse_phase(k)*(H_sup_Ip+H_sup_Im)/sqrt(2) + sin_pulse_phase(k)*(H_sup_Ip-H_sup_Im)/1i/sqrt(2)) + ...  %pulse on carbon
                        2*pi*deltaf*H_sup_offres + 2*pi*methinOffset*H_sup_Sz_CH + ...                         %frequency offset
                        2*pi*H_sup_Jcoupling;                                                                %coupling
                    if (j+3 == 4 && deltaf == min(fRange))
                        M_C = perm_mat_to_make_block_diag(H_sup);
                    end
                    block_C = gpuArray(M_C*(-1i*H_sup - R)*M_C');
                    block1_C = block_C(1:64,1:64);
                    block2_C = block_C(65:end,65:end);
                    block1_C_exp = (expm((block1_C*tau)));
                    block2_C_exp = (expm((block2_C*tau)));
                    prop_C = blkdiag(block1_C_exp,block2_C_exp);
                    rho = (M_C'*prop_C*M_C)*rho;                    
                    % free-precession during transmit channel change
                    rho = prop_free*rho;                                       
                    k = k + 1;
                end
                signal_tmp = gather(rho(2)+1i*rho(3));
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
    
    fitness(i)	= 1e3*(-abs(signal))
end
toc
end












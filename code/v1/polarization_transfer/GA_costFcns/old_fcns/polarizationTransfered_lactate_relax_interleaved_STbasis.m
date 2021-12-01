function[fitness] = polarizationTransfered_lactate_relax_interleaved_STbasis(population)
tic

%% load pre-generated Hamiltonian superoperators for the faseter computation
% only those components are needed that are present in the
% Hamiltonian...taking advantage of linearity
genH_sup = false;
construct_ISTO = false;
if ~genH_sup
    path = 'C:\Users\Somai01\Documents\NMR_stuff\INEPT_optimization\matlabDATAfiles';
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

[npop,nvar]     = size(population);
fitness         = zeros(npop,1);
gamma_C         = 10.71e6;                 % 13C gyromagnetic ratio
gamma_H         = 42.57e6;                 % 1H gyromagnetic ratio
omega_H         = 300e6;                   % 1H resonance frequency for offresonance calculation
omega_C         = omega_H/gamma_H*gamma_C; % 13C resonance frequency for offresonance calculation
J_CH3           = 4.1;                     % coupling in Hz
J_CH            = 3.1;                     % coupling in Hz
J_HH3           = 7;                       % coupling in Hz
B1_range        = 0.5:0.25:1.25;            % pulse amplitude range for compensation
fRange          = -30:15:30;               % offresonance range for compensation
channelChange   = 100e-6;                    % time required to switch between channels on the PIN diode
methinOffset    = 2.6*300;                 % offset of the methin proton relative to the methyle protons

% relaxation parameters
T1_I            = 50;                    % 13C longitudinal
T2_I            = 0.3;                     % 13C transverse
T1_S            = 1.73;                     % 1H longitudinal
T2_S            = 0.1;                     % 1H transverse
% relaxation rates
R1_I = 1/T1_I; R2_I = 1/T2_I; R1_S = 1/T1_S; R2_S = 1/T2_S;

% **********************************************************
% ****** current desing works for I1S3+S spin systems ******
% **********************************************************
n_of_I_spins = 1;
n_of_S_spins = 4;
nspins = n_of_I_spins + n_of_S_spins;

% relaxation rates stored in a matrix for the product operator rate
% calculations
relaxation_rates = repmat([0 R2_I R1_I R2_I],[n_of_I_spins,1]);
relaxation_rates = [relaxation_rates;repmat([0 R2_S R1_S R2_S],[n_of_S_spins,1])];

% Define Pauli matrices
sigma_x =  ([0 1/2; 1/2 0]);
sigma_y =  ([0 -1i/2; 1i/2 0]);
sigma_z =  ([1/2 0; 0 -1/2]);
unit    =  ([1 0; 0 1]);

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

% Carbon spin operators
Ix = L(:,:,1:n_of_I_spins,2);
Iy = L(:,:,1:n_of_I_spins,3);
Iz = L(:,:,1:n_of_I_spins,4);
% Methyle spin operators
Sx_CH3_1 = L(:,:,n_of_I_spins+1,2); Sx_CH3_2 = L(:,:,n_of_I_spins+2,2); Sx_CH3_3 = L(:,:,n_of_I_spins+3,2);
Sy_CH3_1 = L(:,:,n_of_I_spins+1,3); Sy_CH3_2 = L(:,:,n_of_I_spins+2,3); Sy_CH3_3 = L(:,:,n_of_I_spins+3,3);
Sz_CH3_1 = L(:,:,n_of_I_spins+1,4); Sz_CH3_2 = L(:,:,n_of_I_spins+2,4); Sz_CH3_3 = L(:,:,n_of_I_spins+3,4);
% Methin spin operators
Sx_CH = L(:,:,nspins,2);
Sy_CH = L(:,:,nspins,3);
Sz_CH = L(:,:,nspins,4);

% Use magnetic equivalence
Sx_CH3 = Sx_CH3_1+Sx_CH3_2+Sx_CH3_3;
Sy_CH3 = Sy_CH3_1+Sy_CH3_2+Sy_CH3_3;
Sz_CH3 = Sz_CH3_1+Sz_CH3_2+Sz_CH3_3;

% Constuct RL operators
Ip = 1/sqrt(2)*(Ix+1i*Iy);
Im = 1/sqrt(2)*(Ix-1i*Iy);
Sp_CH = -1/sqrt(2)*(Sx_CH+1i*Sy_CH);
Sm_CH = 1/sqrt(2)*(Sx_CH-1i*Sy_CH);
Sp_CH3_1 = -1/sqrt(2)*(Sx_CH3_1+1i*Sy_CH3_1); Sp_CH3_2 = -1/sqrt(2)*(Sx_CH3_2+1i*Sy_CH3_2); Sp_CH3_3 = -1/sqrt(2)*(Sx_CH3_3+1i*Sy_CH3_3);
Sm_CH3_1 = 1/sqrt(2)*(Sx_CH3_1-1i*Sy_CH3_1); Sm_CH3_2 = 1/sqrt(2)*(Sx_CH3_2-1i*Sy_CH3_2); Sm_CH3_3 = 1/sqrt(2)*(Sx_CH3_3-1i*Sy_CH3_3);

Sm_CH3 = Sm_CH3_1+Sm_CH3_2+Sm_CH3_3;
Sp_CH3 = Sp_CH3_1+Sp_CH3_2+Sp_CH3_3;


if construct_ISTO
    %% Build Spherical Tensor operator basis

    % do the housekeeping with L array, i.e. change Ix,Iy etc to Ip,Im etc...
    L(:,:,1:n_of_I_spins,2) = repmat(Im,[1,1,n_of_I_spins,1]);
    L(:,:,1:n_of_I_spins,3) = repmat(Iz,[1,1,n_of_I_spins,1]);
    L(:,:,1:n_of_I_spins,4) = repmat(Ip,[1,1,n_of_I_spins,1]);
    L(:,:,n_of_I_spins+1,2) = Sm_CH3_1;
    L(:,:,n_of_I_spins+1,3) = Sz_CH3_1;
    L(:,:,n_of_I_spins+1,4) = Sp_CH3_1;
    L(:,:,n_of_I_spins+2,2) = Sm_CH3_2;
    L(:,:,n_of_I_spins+2,3) = Sz_CH3_2;
    L(:,:,n_of_I_spins+2,4) = Sp_CH3_2;
    L(:,:,n_of_I_spins+3,2) = Sm_CH3_3;
    L(:,:,n_of_I_spins+3,3) = Sz_CH3_3;
    L(:,:,n_of_I_spins+3,4) = Sp_CH3_3;
    L(:,:,nspins,2) = Sm_CH;
    L(:,:,nspins,3) = Sz_CH;
    L(:,:,nspins,4) = Sp_CH;
    
    T = cell(nspins+1,2*nspins+1,1,1);
    R = zeros(nspins+1,2*nspins+1,1,1);
    for cc1 = 1:nspins+1
        for cc2 = 1:2*nspins+1
            T(cc1,cc2,1) = {zeros(2^nspins,2^nspins)};
        end
    end
    
    T(1,1,1,1,1) = {eye(2^nspins)};
    R(1,1,1,1,1) = relaxation_rates(1,1);
    T(2,1:3,1,1,1) = {Im,Iz,Ip};
    R(2,1:3,1,1,1) = [relaxation_rates(1,2),relaxation_rates(1,3),relaxation_rates(1,4)];
    for n = 2:nspins
        N_STO = sum(2*([0:n-1])+1);                 % number of already generated STOs when adding the nth spin
        % zero out T_tmp entries
        %         for cc1 = 1:nspins+1
        %             for cc2 = 1:2*nspins+1
        %                 for cc3 = 1:nspins+1
        %                     for cc4 = 1:2
        %                         for cc5 = 1:size(cell2mat(T),3)
        %                             T_tmp(cc1,cc2,cc3,cc4,cc5) = {zeros(2^nspins,2^nspins)};
        %                         end
        %                     end
        %                 end
        %             end
        %         end
        T_tmp = cell(nspins+1,2*nspins+1,nspins+1,2,size(cell2mat(T),3));
        R_tmp = zeros(nspins+1,2*nspins+1,nspins+1,2,size(cell2mat(T),3));
        index_accessed = zeros(nspins+1,2*nspins+1,nspins+1,2,size(cell2mat(T),3));
        T_tmp = reshape(T_tmp,[size(T_tmp,1)*size(T_tmp,2)*size(T_tmp,3)*size(T_tmp,4)*size(T_tmp,5),1]);
        for  ttt = 1:size(T_tmp,1)
            T_tmp(ttt) = {zeros(2^nspins,2^nspins)};
        end
        T_tmp = reshape(T_tmp,[nspins+1,2*nspins+1,nspins+1,2,size(T,3)]);
        % loop through all the already generated operators
        for i = 1:N_STO
            Nredundancy_r = size(cell2mat(T),3);
            for r = 1:Nredundancy_r
                q = sum(cumsum(2*([0:n-1])+1) < i);     % gather spin tensor rank
                k = i - sum(2*([0:q-1])+1) - (q+1);     % gather z component
                % k values for a single spin are -1,0,1 therefore 3 terms in the sum
                
                % The added spin operator is E (unity)
                % rank is unchanged
                index_accessed(q+1,k+q+1,q+1,1,r) = index_accessed(q+1,k+q+1,q+1,1,r) + 1;
                idx_tmp = index_accessed(q+1,k+q+1,q+1,1,r);
                T_tmp(q+1,k+q+1,q+1,1,r) = {cell2mat(T_tmp(q+1,k+q+1,q+1,1,r)) + cell2mat(T(q+1,k+q+1,r))*L(:,:,n,1)*1};
                R_tmp(q+1,k+q+1,q+1,1,r) = (idx_tmp-1)/idx_tmp*R_tmp(q+1,k+q+1,q+1,1,r) + 1/idx_tmp*(R(q+1,k+q+1,r)+relaxation_rates(1,1));
                % The added spin operator is Im
                if q > 0
                    % rank is decreased
                    if q - 1 >= abs(k-1) % check whether the decreased rank tensor is allowed (because -k=<q<=k)
                        index_accessed(q-1+1,k-1+q-1+1,q+1,2,r) = index_accessed(q-1+1,k-1+q-1+1,q+1,2,r) + 1;
                        idx_tmp = index_accessed(q-1+1,k-1+q-1+1,q+1,2,r);
                        T_tmp(q-1+1,k-1+q-1+1,q+1,2,r) = {cell2mat(T_tmp(q-1+1,k-1+q-1+1,q+1,2,r)) + cell2mat(T(q+1,k+q+1,r))*L(:,:,n,2)*CG_coeff(q,k,1,-1,q-1,k-1)};
                        R_tmp(q-1+1,k-1+q-1+1,q+1,2,r) = (idx_tmp-1)/idx_tmp*R_tmp(q-1+1,k-1+q-1+1,q+1,2,r) + 1/idx_tmp*(R(q+1,k+q+1,r)+relaxation_rates(n,2))*CG_coeff(q,k,1,-1,q-1,k-1);
                    end
                    % rank is unchanged
                    if q >= abs(k-1)% check whether the decreased rank tensor is allowed (because -k=<q<=k)
                        index_accessed(q+1,k-1+q+1,q+1,2,r) = index_accessed(q+1,k-1+q+1,q+1,2,r) + 1;
                        idx_tmp = index_accessed(q+1,k-1+q+1,q+1,2,r);
                        T_tmp(q+1,k-1+q+1,q+1,2,r) = {cell2mat(T_tmp(q+1,k-1+q+1,q+1,2,r)) + cell2mat(T(q+1,k+q+1,r))*L(:,:,n,2)*CG_coeff(q,k,1,-1,q,k-1)};
                        R_tmp(q+1,k-1+q+1,q+1,2,r) = (idx_tmp-1)/idx_tmp*R_tmp(q+1,k-1+q+1,q+1,2,r) +1/idx_tmp*(R(q+1,k+q+1,r))+relaxation_rates(n,2)*CG_coeff(q,k,1,-1,q,k-1);
                    end
                    % rank is increased
                    index_accessed(q+1+1,k-1+q+1+1,q+1,2,r) = index_accessed(q+1+1,k-1+q+1+1,q+1,2,r) + 1;
                    idx_tmp = index_accessed(q+1+1,k-1+q+1+1,q+1,2,r);
                    T_tmp(q+1+1,k-1+q+1+1,q+1,2,r) = {cell2mat(T_tmp(q+1+1,k-1+q+1+1,q+1,2,r)) + cell2mat(T(q+1,k+q+1,r))*L(:,:,n,2)*CG_coeff(q,k,1,-1,q+1,k-1)};
                    R_tmp(q+1+1,k-1+q+1+1,q+1,2,r) = (idx_tmp-1)/idx_tmp*R_tmp(q+1+1,k-1+q+1+1,q+1,2,r) + 1/idx_tmp*(R(q+1,k+q+1,r))+relaxation_rates(n,2)*CG_coeff(q,k,1,-1,q+1,k-1);
                else% if the rank was zero, it only can increase (triangle ineq)
                    index_accessed(q+1+1,k-1+q+1+1,q+1,2,r) = index_accessed(q+1+1,k-1+q+1+1,q+1,2,r) + 1;
                    idx_tmp = index_accessed(q+1+1,k-1+q+1+1,q+1,2,r);
                    T_tmp(q+1+1,k-1+q+1+1,q+1,2,r) = {cell2mat(T_tmp(q+1+1,k-1+q+1+1,q+1,2,r)) + cell2mat(T(q+1,k+q+1,r))*L(:,:,n,2)*CG_coeff(q,k,1,-1,q+1,k-1)};
                    R_tmp(q+1+1,k-1+q+1+1,q+1,2,r) = (idx_tmp-1)/idx_tmp*R_tmp(q+1+1,k-1+q+1+1,q+1,2,r) + 1/idx_tmp*(R(q+1,k+q+1,r))+relaxation_rates(n,2)*CG_coeff(q,k,1,-1,q+1,k-1);
                end
                % The added spin operator is Iz
                if q > 0
                    % rank is decreased
                    if q - 1 >= abs(k) % check whether the decreased rank tensor is allowed (because -k=<q<=k)
                        index_accessed(q-1+1,k+q-1+1,q+1,2,r) = index_accessed(q-1+1,k+q-1+1,q+1,2,r) + 1;
                        idx_tmp = index_accessed(q-1+1,k+q-1+1,q+1,2,r);
                        T_tmp(q-1+1,k+q-1+1,q+1,2,r) = {cell2mat(T_tmp(q-1+1,k+q-1+1,q+1,2,r)) + cell2mat(T(q+1,k+q+1,r))*L(:,:,n,3)*CG_coeff(q,k,1,0,q-1,k)};
                        R_tmp(q-1+1,k+q-1+1,q+1,2,r) = (idx_tmp-1)/idx_tmp*R_tmp(q-1+1,k+q-1+1,q+1,2,r) + 1/idx_tmp*(R(q+1,k+q+1,r))+relaxation_rates(n,3)*CG_coeff(q,k,1,0,q-1,k);
                    end
                    % rank is unchanged
                    index_accessed(q+1,k+q+1,q+1,2,r) = index_accessed(q+1,k+q+1,q+1,2,r) + 1;
                    idx_tmp = index_accessed(q+1,k+q+1,q+1,2,r);
                    T_tmp(q+1,k+q+1,q+1,2,r) = {cell2mat(T_tmp(q+1,k+q+1,q+1,2,r)) + cell2mat(T(q+1,k+q+1,r))*L(:,:,n,3)*CG_coeff(q,k,1,0,q,k)};
                    R_tmp(q+1,k+q+1,q+1,2,r) = (idx_tmp-1)/idx_tmp*R_tmp(q+1,k+q+1,q+1,2,r) + 1/idx_tmp*(R(q+1,k+q+1,r))+relaxation_rates(n,3)*CG_coeff(q,k,1,0,q,k);
                    % rank is increased
                    index_accessed(q+1+1,k+q+1+1,q+1,2,r) = index_accessed(q+1+1,k+q+1+1,q+1,2,r) + 1;
                    idx_tmp = index_accessed(q+1+1,k+q+1+1,q+1,2,r);
                    T_tmp(q+1+1,k+q+1+1,q+1,2,r) = {cell2mat(T_tmp(q+1+1,k+q+1+1,q+1,2,r)) + cell2mat(T(q+1,k+q+1,r))*L(:,:,n,3)*CG_coeff(q,k,1,0,q+1,k)};
                    R_tmp(q+1+1,k+q+1+1,q+1,2,r) = (idx_tmp-1)/idx_tmp*R_tmp(q+1+1,k+q+1+1,q+1,2,r) + 1/idx_tmp*(R(q+1,k+q+1,r))+relaxation_rates(n,3)*CG_coeff(q,k,1,0,q+1,k);
                else% if the rank was zero, it only can increase (triangle ineq)
                    index_accessed(q+1+1,k+q+1+1,q+1,2,r) = index_accessed(q+1+1,k+q+1+1,q+1,2,r) + 1;
                    idx_tmp = index_accessed(q+1+1,k+q+1+1,q+1,2,r);
                    T_tmp(q+1+1,k+q+1+1,q+1,2,r) = {cell2mat(T_tmp(q+1+1,k+q+1+1,q+1,2,r)) + cell2mat(T(q+1,k+q+1,r))*L(:,:,n,3)*CG_coeff(q,k,1,0,q+1,k)};
                    R_tmp(q+1+1,k+q+1+1,q+1,2,r) = (idx_tmp-1)/idx_tmp*R_tmp(q+1+1,k+q+1+1,q+1,2,r) + 1/idx_tmp*(R(q+1,k+q+1,r))+relaxation_rates(n,3)*CG_coeff(q,k,1,0,q+1,k);
                end
                % The added spin operator is Ip
                if q > 0
                    % rank is decreased
                    if q - 1 >= abs(k+1) % check whether the decreased rank tensor is allowed (because -k=<q<=k)
                        index_accessed(q-1+1,k+1+q-1+1,q+1,2,r) = index_accessed(q-1+1,k+1+q-1+1,q+1,2,r) + 1;
                        idx_tmp = index_accessed(q-1+1,k+1+q-1+1,q+1,2,r);
                        T_tmp(q-1+1,k+1+q-1+1,q+1,2,r) = {cell2mat(T_tmp(q-1+1,k+1+q-1+1,q+1,2,r)) + cell2mat(T(q+1,k+q+1,r))*L(:,:,n,4)*CG_coeff(q,k,1,1,q-1,k+1)};
                        R_tmp(q-1+1,k+1+q-1+1,q+1,2,r) = (idx_tmp-1)/idx_tmp*R_tmp(q-1+1,k+1+q-1+1,q+1,2,r) + 1/idx_tmp*(R(q+1,k+q+1,r))+relaxation_rates(n,4)*CG_coeff(q,k,1,1,q-1,k+1);
                    end
                    % rank is unchanged
                    if q >= abs(k+1) % check whether the decreased rank tensor is allowed (because -k=<q<=k)
                        index_accessed(q+1,k+1+q+1,q+1,2,r) = index_accessed(q+1,k+1+q+1,q+1,2,r) + 1;
                        idx_tmp = index_accessed(q+1,k+1+q+1,q+1,2,r);
                        T_tmp(q+1,k+1+q+1,q+1,2,r) = {cell2mat(T_tmp(q+1,k+1+q+1,q+1,2,r)) + cell2mat(T(q+1,k+q+1,r))*L(:,:,n,4)*CG_coeff(q,k,1,1,q,k+1)};
                        R_tmp(q+1,k+1+q+1,q+1,2,r) = (idx_tmp-1)/idx_tmp*R_tmp(q+1,k+1+q+1,q+1,2,r) + 1/idx_tmp*(R(q+1,k+q+1,r))+relaxation_rates(n,4)*CG_coeff(q,k,1,1,q,k+1);
                    end
                    % rank is increased
                    index_accessed(q+1+1,k+1+q+1+1,q+1,2,r) = index_accessed(q+1+1,k+1+q+1+1,q+1,2,r) + 1;
                    idx_tmp = index_accessed(q+1+1,k+1+q+1+1,q+1,2,r);
                    T_tmp(q+1+1,k+1+q+1+1,q+1,2,r) = {cell2mat(T_tmp(q+1+1,k+1+q+1+1,q+1,2,r)) + cell2mat(T(q+1,k+q+1,r))*L(:,:,n,4)*CG_coeff(q,k,1,1,q+1,k+1)};
                    R_tmp(q+1+1,k+1+q+1+1,q+1,2,r) = (idx_tmp-1)/idx_tmp*R_tmp(q+1+1,k+1+q+1+1,q+1,2,r) + 1/idx_tmp*(R(q+1,k+q+1,r))+relaxation_rates(n,4)*CG_coeff(q,k,1,1,q+1,k+1);
                else% if the rank was zero, it only can increase (triangle ineq)
                    index_accessed(q+1+1,k+1+q+1+1,q+1,2,r) = index_accessed(q+1+1,k+1+q+1+1,q+1,2,r) + 1;
                    idx_tmp = index_accessed(q+1+1,k+1+q+1+1,q+1,2,r);
                    T_tmp(q+1+1,k+1+q+1+1,q+1,2,r) = {cell2mat(T_tmp(q+1+1,k+1+q+1+1,q+1,2,r)) + cell2mat(T(q+1,k+q+1,r))*L(:,:,n,4)*CG_coeff(q,k,1,1,q+1,k+1)};
                    R_tmp(q+1+1,k+1+q+1+1,q+1,2,r) = (idx_tmp-1)/idx_tmp*R_tmp(q+1+1,k+1+q+1+1,q+1,2,r) + 1/idx_tmp*(R(q+1,k+q+1,r))+relaxation_rates(n,4)*CG_coeff(q,k,1,1,q+1,k+1);
                end
            end
            fprintf("spin %i rank %i component %i done\n",n,q,k)
        end
        T = T_tmp;
        R = R_tmp;
        T = reshape(T,[nspins+1,2*nspins+1,size(T,3)*size(T,4)*size(T,5)]);
        R = reshape(R,[nspins+1,2*nspins+1,size(R,3)*size(R,4)*size(R,5)]);
    end
    
    B = cell(size(T,1)*size(T,2)*size(T,3),1);
    R_sup = zeros(size(T,1)*size(T,2)*size(T,3),1);
    basis_index = 1;
    for cc1 = 1:nspins+1
        for cc2 = 1:2*nspins+1
            for cc3 = 1:size(T,3)
                if (sum(sum(cell2mat(T(cc1,cc2,cc3)) ~= zeros(2^nspins,2^nspins))) > 0)
                    B((cc1-1)*(2*nspins+1)*size(T,3) + (cc2-1)*size(T,3) + cc3) = T(cc1,cc2,cc3);
                    R_sup((cc1-1)*(2*nspins+1)*size(T,3) + (cc2-1)*size(T,3) + cc3) = R(cc1,cc2,cc3);
                end
            end
        end
    end
    indeces2keep = ~cellfun('isempty',B);
    B = B(indeces2keep);
    R_sup = R_sup(indeces2keep);
    disp('Tensor generation done')
    for ttt = 1:size(B,1)
        for ttt2 = ttt+1:size(B,1)
            if (sum(sum(cell2mat(B(ttt)) == cell2mat(B(ttt2)))) == (2^nspins*2^nspins))
                B(ttt2) = {zeros(2^nspins,2^nspins)};
                R_sup(ttt2) = 0;
            end
        end
    end
    disp('Redundancy test done')
    for ttt = 1:size(B,1)
        if (sum(sum(cell2mat(B(ttt)) == zeros(2^nspins,2^nspins))) == (2^nspins*2^nspins))
            B(ttt) = {[]};
        end
    end
    B = B(~cellfun('isempty',B));
    R_sup = R_sup(~cellfun('isempty',B));
    disp('Redundant entries dropped')
    
    save(fullfile('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization',strcat('ISTO_',num2str(nspins),'spins')),'B')
    save(fullfile('C:\Users\Somai01\Documents\pulseSequenceDesign\NMR_stuff\INEPT_optimization',strcat('R_sup_',num2str(nspins),'spins')),'R_sup')
else
    load(fullfile(path,strcat('ISTO_',num2str(nspins),'spins')))
    load(fullfile(path,strcat('R_sup_',num2str(nspins),'spins')))
    R_sup = diag(R_sup);
end

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
    return;
end

H_sup_Sm = (H_sup_Sm_CH+H_sup_Sm_CH3);
H_sup_Sp = (H_sup_Sp_CH+H_sup_Sp_CH3);
H_sup_Sz = (H_sup_Sz_CH+H_sup_Sz_CH3);
H_sup_Jcoupling = (J_CH3*H_sup_IzSz_CH3+J_CH*H_sup_IzSz_CH+J_HH3*H_sup_Sz_CHSz_CH3);
H_sup_offres = (H_sup_Iz + gamma_H/gamma_C*(H_sup_Sz_CH+H_sup_Sz_CH3));

%% fitness evaluation
for i = npop:-1:1
    parameters = population(i,1:nvar);
    signal = zeros(length(B1_range),length(B1_range));
    contamination = 0;
    cos_pulse_phase = cos(parameters((1:3:nvar)+1));
    sin_pulse_phase = sin(parameters((1:3:nvar)+1));
    signal_tmp = 0;
    for deltaf = fRange
        H_sup_free = 2*pi*deltaf*H_sup_offres + 2*pi*methinOffset*H_sup_Sz_CH + ...                %frequency offset
            2*pi*H_sup_Jcoupling;                                                                %coupling
        prop_free = (expm(gpuArray((-1i*H_sup_free - R_sup)*channelChange)));
        bb_C = 1;
        sum_weight = 0;
        for B1_C = B1_range
            bb_H = 1;
            for B1_H = B1_range
                % Initial state = spins along z-axis, polariation on protons
                rho = zeros((2^nspins)^2,1,'gpuArray');
%                 rho(7) = 1;
%                 rho(10) = 1;
%                 rho(13) = 1;
%                 rho(16) = 1;
                rho(4) = 1;
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
%                     if (j == 1 && deltaf == min(fRange))
%                         M_H = perm_mat_to_make_block_diag(H_sup);
%                         R1 = R_sup(1:512,1:512);
%                         R2 = R_sup(513:end,513:end);
%                     end
%                     block_H = gpuArray(M_H*(-1i*H_sup - 0*R_sup)*M_H');
%                     block1_H = block_H(1:512,1:512);
%                     block2_H = block_H(513:end,513:end);
%                     block1_H_exp = (expm((block1_H*tau)));
%                     block2_H_exp = (expm((block2_H*tau)));
%                     prop_H = blkdiag(block1_H_exp,block2_H_exp);
%                     rho = (M_H'*prop_H*M_H)*rho;                                                       %build propagator
                   %prop_H = expm(gpuArray(-1i*H_sup - R_sup)*tau);
                    %rho = prop_H*rho;
                    rho = expv(tau,gather(-1i*H_sup - R_sup),gather(rho),1e-7,30);

                    
                    % free-precession during transmit channel change
                    rho = prop_free*rho;
                    
                    % carbon
                    %j = j + 3;
                    k = k + 1;
                    B1 = B1_C*parameters(j+3);
                    tau = parameters(j+2+3);
                    H_sup = 2*pi*gamma_C*B1*(cos_pulse_phase(k)*(H_sup_Ip+H_sup_Im)/sqrt(2) + sin_pulse_phase(k)*(H_sup_Ip-H_sup_Im)/1i/sqrt(2)) + ... %pulse on carbon
                        2*pi*deltaf*H_sup_offres + 2*pi*methinOffset*H_sup_Sz_CH + ...                         %frequency offset
                        2*pi*H_sup_Jcoupling;                                                                %coupling
%                     if (j == 1 && deltaf == min(fRange))
%                         M_C = perm_mat_to_make_block_diag(H_sup);
%                     end
%                     block_C = gpuArray(M_C*(-1i*H_sup - 0*R_sup)*M_C');
%                     block1_C = block_C(1:64,1:64);
%                     block2_C = block_C(65:end,65:end);
%                     block1_C_exp = (expm((block1_C*tau)));
%                     block2_C_exp = (expm((block2_C*tau)));
%                     prop_C = blkdiag(block1_C_exp,block2_C_exp);
%                     rho = (M_C'*prop_C*M_C)*rho;
                      %prop_C = expm(gpuArray(-1i*H_sup - R_sup)*tau);
                      %rho = prop_C*rho;
                      rho = expv(tau,gather(-1i*H_sup - R_sup),gather(rho),1e-7,30);
                    
                    % free-precession during transmit channel change
                    rho = prop_free*rho;
                    k = k + 1;
                end
                %signal_tmp = gather(rho(2)+1i*rho(3));
                signal_tmp = gather(rho(5)+1i*rho(6)+rho(8)+1i*rho(9)+rho(11)+1i*rho(12));
                weight = 1;%/(1+(1-B1_C)^2+(1-B1_H)^2);
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












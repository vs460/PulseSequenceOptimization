clear all; close all
%% load pulse sequence
[file,path] = uigetfile({'*.*'});
load(fullfile(path,file))
%% initialize parameters
parameters    = pars.parameters;
RF            = pars.RF;
nvar          = length(parameters)/3;
gamma_C       = 10.71e6;                 % 13C gyromagnetic ratio
gamma_H       = 42.57e6;                 % 1H gyromagnetic ratio
omega_H       = 299.709e6;               % 1H resonance frequency for offresonance calculation
omega_C       = omega_H/gamma_H*gamma_C; % 13C resonance frequency for offresonance calculation
J_CH3         = 4.1;                     % coupling in Hz
J_CH          = 146;                    % coupling in Hz
J_HH3         = 7.1;                     % coupling in Hz
B1_range_X    = linspace(1,1,1);%pars.B1_range;           % pulse amplitude range for compensation
fRange        = linspace(-0,0,1);%pars.fRange;             % offresonance range for compensation
metOffset     = -2.79*299.709;
min_C_amp     = 0;
min_H_amp     = 0;
H_maxAmp      = 0.25e-4;
C_maxAmp      = 1e-4;
lB1r_X        = length(B1_range_X);
lB0r          = length(fRange);
channelChange = 200e-6;

% B0 inhom distribution
LW = max(fRange)-min(fRange);
lineBroadening = normpdf(fRange,0,LW/2.355);
lineBroadening = lineBroadening/max(lineBroadening);
%B1weight_X = normpdf(log(B1_range_X)-0.2,0,4/2.355);
B1weight_X = normpdf(B1_range_X-1,0,1.5/2.355);
B1weight_X = B1weight_X/max(B1weight_X);

% relaxation parameters
T1_I = 7.2;             % 13C longitudinal  [s]
T2_I = 0.3;             % 13C transverse    [s]
T1_S = 1.73;            % 1H longitudinal   [s]
T2_S = 0.1;             % 1H transverse     [s]
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
sigma_x = ([0 1/2; 1/2 0]);
sigma_y = ([0 -1i/2; 1i/2 0]);
sigma_z = ([1/2 0; 0 -1/2]);
unit    = ([1 0; 0 1]);

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

Sx_CH3 = sum(L(:,:,n_of_I_spins+1:nspins-1,2),3);
Sy_CH3 = sum(L(:,:,n_of_I_spins+1:nspins-1,3),3);
Sz_CH3 = sum(L(:,:,n_of_I_spins+1:nspins-1,4),3);

Sx_CH = L(:,:,nspins,2);
Sy_CH = L(:,:,nspins,3);
Sz_CH = L(:,:,nspins,4);


% construct cartesian product operator basis
% basis elements are represented as a 1x4 row vector:
% entry 1 = E, 2 = Lx, 3 = Ly, 4 = Lz
% and the basis operators are the product of the permutations (with repetitions) of these
for bb = 1:nspins
    permutation_mtx_tmp(:,bb) = reshape(repmat(mod((1:4^(nspins-bb+1))-1,4)'+1,[1,4^(bb-1)])',[1,4^nspins])';
end
permutation_mtx = permutation_mtx_tmp;

% constructing cell with the basis operators and the relaxation supermatrix
B             = cell((2^nspins)^2,1);
R             = zeros((2^nspins)^2,1);
RelMtx        = zeros(size(Ix));
RelMtx_struct = zeros(size(Ix));

B_block = [];
for bb = 1:4^nspins
    basis_prod_tmp = eye(2^nspins);
    for ns = 1:nspins
        basis_prod_tmp = basis_prod_tmp*L(:,:,ns,permutation_mtx(bb,ns)); % basis operator is the product of 5 terms according to the correspopnding permutation
        R(bb) = R(bb) + relaxation_rates(ns,permutation_mtx(bb,ns));      % relaxation rate of product operators = sum of correspoding single rates
    end
    if bb >16    % normalization according to spindynamica exported basis set
        B(bb) = {basis_prod_tmp*(2^(sum(permutation_mtx(bb,:)>1)-2))/sqrt(2)};
    elseif bb > 1
        B(bb) = {basis_prod_tmp/2/sqrt(2)};
    else
        B(bb) = {basis_prod_tmp/4/sqrt(2)};
    end
    RelMtx = RelMtx + R(bb)*(cell2mat(B(bb))~=0);
    RelMtx_struct = RelMtx_struct + (cell2mat(B(bb))~=0);
end
RelMtx = (RelMtx./RelMtx_struct);

Sx          = (Sx_CH+Sx_CH3);
Sy          = (Sy_CH+Sy_CH3);
Sz          = (Sz_CH+Sz_CH3);
IzSz_CH3    = (Iz*Sz_CH3);
IzSz_CH     = (Iz*Sz_CH);
Sz_CHSz_CH3 = (Sz_CH*Sz_CH3);
Jcoupling   = (J_CH3*IzSz_CH3+J_CH*IzSz_CH+J_HH3*(Sz_CHSz_CH3+0*Sy_CH*Sy_CH3+0*Sx_CH*Sx_CH3));
offres      = (Iz + gamma_H/gamma_C*(Sz_CH+Sz_CH3));
Idmtx       = eye(2^nspins);

%% initialization
initState   = Iz;
traceNorm   = trace(Iz'*Iz);
coil        = Ix+1i*Iy;
%% fitness evaluation
sum_weight  = 0;
signal      = zeros(lB1r_X,lB0r);
signal_Z    = zeros(lB1r_X,lB0r);
rho_array   = zeros(lB1r_X,lB0r,(2^nspins)^2); 
h = waitbar(0,'Evaluation over grid');
for ff = 1:length(fRange)
    Wf      = lineBroadening(ff);
    df      = fRange(ff);
    H_free  = 2*pi*df*offres + 2*pi*metOffset*Sz_CH3 + ...                %frequency offset
        2*pi*Jcoupling;
    parfor bb_C = 1:length(B1_range_X)
        B1_C = B1_range_X(bb_C)
        WBC  = B1weight_X(bb_C);
        rho  = initState(:);
        %%%%%%%%%% Pulse sequence %%%%%%%%%%
        for j = 1:nvar
            B1_I = B1_C*RF.Amp(j);
            ph   = RF.Phase(j);
            tau  = RF.Taus(j);
            
            %%% 1H-Xnuc pulse %%%
            H   = 2*pi*gamma_C*B1_I*(cos(ph)*Ix + sin(ph)*Iy)+...
                H_free;
            PH  = kron(eye(2^nspins),H)-kron(H.',eye(2^nspins));
            L   = -1i*PH - diag(RelMtx(:));
            rho = expv(tau,L,rho,1e-7,30);
            
            
        end
        rho_array(bb_C,ff,:) = rho;
        rho                  = reshape(rho,[size(initState)]);
        weight               = WBC*Wf;
        signal_tmp           = trace(coil'*rho)/traceNorm;
        signal(bb_C,ff) = signal_tmp;
        signal_Z(bb_C,ff) = trace(Sz_CH'*rho)/traceNorm;;
        sum_weight           = sum_weight + weight;
    end
    waitbar(((ff-1)*lB1r_X+lB1r_X)/(lB0r*lB1r_X))
    
end
pars.signal = signal;
signal_tred  = sum(sum(sum(signal)))/sum_weight;
avgTransfer  = signal_tred

%% detection
[~,I1] = min(abs(B1_range_X - 1));
[~,I2] = min(abs(fRange-0));
rho          = reshape(rho_array(I1,I2,:),[size(initState)]);
peakTransfer = trace(coil'*rho)/traceNorm
nsteps       = 4096;
BW           = 4000;
f            = linspace(-BW/2,BW/2,nsteps);
H_free       = 2*pi*metOffset*Sz_CH3 + ...
    2*pi*Jcoupling;
PH_free      = kron(eye(2^nspins),H_free)-kron(H_free.',eye(2^nspins));
L_free       = -1i*PH_free - diag(RelMtx(:));
prop_free    = expm(L_free*1/BW);

for n = 1:nsteps
    fid(n)= trace(coil'*rho)/traceNorm;
    rho   = prop_free*rho(:);
    rho   = reshape(rho,[2^nspins,2^nspins]);
end
%%
fid_apd   = fid.*exp(-0.005*linspace(1,nsteps,nsteps));
spectrum  = ifftshift(ifft(fid_apd));
plot(f,abs(spectrum))
sumIntensity = abs(trapz(spectrum))


















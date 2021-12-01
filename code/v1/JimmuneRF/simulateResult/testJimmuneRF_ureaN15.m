%function [] = testJimmuneRF_ureaN15(pars)
 clear all; close all 
%% load pulse sequence (pars.mat file)
[file,path] = uigetfile({'*.*'});
load(fullfile(path,file))
%% initialize parameters
parameters    = pars.parameters;
RF            = pars.RF;
nvar          = length(parameters)/3;
gamma_N       = pars.gX;                 % Xnuc gyromagnetic ratio
gamma_H       = pars.gH;                 % 1H gyromagnetic ratio
J_NH          = -90;                     % coupling in Hz
B1_range_X    = linspace(0.5,2,16);      % pulse amplitude range for compensation
fRange        = linspace(-50,50,100);    % offresonance range for compensation
min_N_amp     = pars.min_X_amp;
N_maxAmp      = pars.X_maxAmp;
lB1r_X        = length(B1_range_X);
lB0r          = length(fRange);
% B0 inhom distribution
LW = (max(fRange)-min(fRange))/2;               % FWHM is assumed to be the half of the frequency range
lineBroadening = normpdf(fRange,0,LW/2.355);    % set variance to match the expected field inhom
lineBroadening = lineBroadening/max(lineBroadening);
% B1 inhom distribution
B1weight_X = normpdf(B1_range_X-1,0,1.5/2.355); % set variance to match the expected coil Tx profile
B1weight_X = B1weight_X/max(B1weight_X);

% relaxation parameters
T1_I = pars.T1_I;         % Xnuc longitudinal
T2_I = pars.T2_I;         % Xnuc transverse
T1_S = pars.T1_S;         % 1H longitudinal
T2_S = pars.T2_S;         % 1H transverse
% relaxation rates
R1_I = 1/T1_I; R2_I = 1/T2_I; R1_S = 1/T1_S; R2_S = 1/T2_S;

% **********************************************************
% ******* current desing works for I1S3 spin systems *******
% **********************************************************
n_of_I_spins = 1;
n_of_S_spins = 2;
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

Sx = sum(L(:,:,n_of_I_spins+1:nspins,2),3);
Sy = sum(L(:,:,n_of_I_spins+1:nspins,3),3);
Sz = sum(L(:,:,n_of_I_spins+1:nspins,4),3);

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
        basis_prod_tmp = basis_prod_tmp*L(:,:,ns,permutation_mtx(bb,ns)); % basis operator is the product of nspins terms
        R(bb) = R(bb) + relaxation_rates(ns,permutation_mtx(bb,ns));      % relaxation rate of product operators = sum of correspoding single rates
    end    % normalization 
    B(bb) = {basis_prod_tmp*(2^(sum(permutation_mtx(bb,:)>1)-2))/sqrt(2)};
    RelMtx = RelMtx + R(bb)*(cell2mat(B(bb))~=0);
    RelMtx_struct = RelMtx_struct + (cell2mat(B(bb))~=0);
end
RelMtx = (RelMtx./RelMtx_struct);

Jcoupling = J_NH*Iz*Sz; % use magnetic equivalence and secular approximation
offres    = (Iz +  gamma_H/gamma_N*Sz);
Idmtx     = eye(2^nspins);
%% initialize states
initState = Iz;
traceNorm = trace(Iz'*Iz); % normalise the result to the initial state
coil      = Ix+1i*Iy;      % detection state
%% fitness evaluation
sum_weight = 0;
signal     = zeros(lB1r_X,lB0r);
signalZ    = zeros(lB1r_X,lB0r);
h = waitbar(0,'Evaluation over grid');
for ff = 1:length(fRange)
    Wf     = lineBroadening(ff);
    df     = fRange(ff);
    % background generator Hamiltonian 
    H_free = 2*pi*df*offres + 2*pi*Jcoupling;
    bb_C   = 1;
    parfor bb_C = 1:length(B1_range_X)
        B1_C = B1_range_X(bb_C)
        WBC  = B1weight_X(bb_C);
        rho  = initState(:);
        %%%%%%%%%% Pulse sequence %%%%%%%%%%
        for j=1:nvar
            B1_I = B1_C*RF.Amp(j);
            ph   = RF.Phase(j);
            tau  = RF.Taus(j);
            
            %%% pulse %%%
            H   = 2*pi*gamma_N*B1_I*(cos(ph)*Ix + sin(ph)*Iy)+...
                H_free;
            PH  = kron(eye(2^nspins),H)-kron(H.',eye(2^nspins));
            L   = -1i*PH - diag(RelMtx(:)); % putting together the Liouvillian
            rho = expv(tau,L,rho,1e-7,30);  % Krylov propagation  
        end
        rho_array(bb_C,ff,:) = rho;
        rho                  = reshape(rho,[size(initState)]);
        weight               = WBC*Wf;
        signal_tmp           = trace(coil'*rho)/traceNorm;
        signal(bb_C,ff)      = signal_tmp;
        signalZ(bb_C,ff)    = trace(Iz'*rho)/traceNorm;
        sum_weight           = sum_weight + weight;
    end
    waitbar(((ff-1)*lB1r_X+lB1r_X)/(lB0r*lB1r_X))    
end
close(h)

pars.signal = signal;
signal_tred  = sum(sum(signal))/sum_weight;
avgTransfer  = signal_tred

%% calculating spectrum corresponding to perfect B0 and B1
[~,I1]       = min(abs(B1_range_X - 1));
[~,I2]       = min(abs(fRange-0));
rho          = reshape(rho_array(I1,I2,:),[size(initState)]);
rhoAvg       = reshape(sum(sum(rho_array,1),2),[size(initState)])/sum_weight;
peakSignal   = abs(trace(coil'*rho)/traceNorm)
nsteps       = 4096;
BW           = 4000;
f            = linspace(-BW/2,BW/2,nsteps);
% free precession propagator
PH_free      = kron(eye(2^nspins),H_free)-kron(H_free.',eye(2^nspins));
L_free       = -1i*PH_free - diag(RelMtx(:));
prop_free    = expm(L_free*1/BW);
% acquiring the FID
for n = 1:nsteps
    fid(n)   = trace(coil'*rho)/traceNorm;
    fidAvg(n)= trace(coil'*rhoAvg)/traceNorm;
    rho      = prop_free*rho(:);
    rhoAvg   = prop_free*rhoAvg(:);
    rho      = reshape(rho,[2^nspins,2^nspins]);
    rhoAvg   = reshape(rhoAvg,[2^nspins,2^nspins]);
end
%%
figure
fid_apd      = fid.*exp(-0.005*linspace(1,nsteps,nsteps));
fidAvg_apd   = fidAvg.*exp(-0.005*linspace(1,nsteps,nsteps));
spectrum     = ifftshift(ifft(fid_apd));
spectrumAvg  = ifftshift(ifft(fidAvg_apd));
plot(f,abs(spectrum))
hold on
plot(f,abs(spectrumAvg))
legend('Perfect spectrum','Ensemble average spectrum')
peakIntegral = abs(trapz(spectrum))











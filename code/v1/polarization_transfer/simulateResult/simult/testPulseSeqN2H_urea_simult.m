clear all; close all
%% load pulse sequence 
[file,path] = uigetfile({'*.*'});
load(fullfile(path,file))
%% initialize parameters
parameters    = pars.parameters;
RF            = pars.RF;
nvar          = length(parameters)/5;
gamma_N       = pars.gX;                 % Xnuc gyromagnetic ratio
gamma_H       = pars.gH;                 % 1H gyromagnetic ratio
omega_H       = 299.709e6;               % 1H resonance frequency for offresonance calculation
omega_N       = omega_H/gamma_H*gamma_N; % 15N resonance frequency for offresonance calculation
J_NH          = -90;                     % coupling in Hz
B1_range_H    = linspace(0.1,3,10);%pars.B1_range;           % pulse amplitude range for compensation
B1_range_X    = linspace(0.1,3,10);%pars.B1_range;           % pulse amplitude range for compensation
fRange        = linspace(-200,200,15);%pars.fRange;             % offresonance range for compensation
min_N_amp     = pars.min_X_amp;
min_H_amp     = pars.min_H_amp;
N_maxAmp      = pars.X_maxAmp;
H_maxAmp      = pars.H_maxAmp;
lB1r_H        = length(B1_range_H);
lB1r_X        = length(B1_range_X);
lB0r          = length(fRange);
presPol       = 0.73;
% B0 inhom distribution
LW = max(fRange)-min(fRange);
lineBroadening = normpdf(fRange,0,LW/2.355);
lineBroadening = lineBroadening/max(lineBroadening);
% B1 inhom distribution
B1weight_H = normpdf(B1_range_H-1,0,1.5/2.355);
B1weight_H = B1weight_H/max(B1weight_H);
%B1weight_X = normpdf(log(B1_range_X)-0.2,0,4/2.355);
B1weight_X = normpdf(B1_range_X-1,0,1.5/2.355);
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
        basis_prod_tmp = basis_prod_tmp*L(:,:,ns,permutation_mtx(bb,ns)); % basis operator is the product of 5 terms according to the correspopnding permutation
        R(bb) = R(bb) + relaxation_rates(ns,permutation_mtx(bb,ns));      % relaxation rate of product operators = sum of correspoding single rates
    end    % normalization according to spindynamica exported basis set
    B(bb) = {basis_prod_tmp*(2^(sum(permutation_mtx(bb,:)>1)-2))/sqrt(2)};
    RelMtx = RelMtx + R(bb)*(cell2mat(B(bb))~=0);
    RelMtx_struct = RelMtx_struct + (cell2mat(B(bb))~=0);
end
RelMtx = (RelMtx./RelMtx_struct);

Jcoupling = J_NH*Iz*Sz; % use magnetic equivalence and secular approximation
offres    = (Iz/(gamma_H/gamma_N) + Sz);
Idmtx     = eye(2^nspins);
%% initialize states
initState     = Iz;
traceNorm     = trace(Iz'*Iz);
coil          = Sx+1i*Sy;
%% fitness evaluation
sum_weight = 0;
signal     = zeros(lB1r_X,lB1r_H,lB0r);
h = waitbar(0,'Evaluation over grid');
for ff = 1:length(fRange)
    Wf     = lineBroadening(ff);
    df     = fRange(ff);
    H_free = 2*pi*df*offres + 2*pi*Jcoupling;
    bb_C   = 1;
    for B1_C = B1_range_X
        bb_H = 1;
        for B1_H = B1_range_H
            rho  = initState(:);
            %%%%%%%%%% Pulse sequence %%%%%%%%%%
            for j = 1:nvar
                B1_S = B1_H*RF.protonAmp(j);
                B1_I = B1_C*RF.XAmp(j);
                %B1_S = min(round(abs(parameters(j))/min_H_amp),1)*sign(parameters(j))*B1_H*max(min(abs(parameters(j)),H_maxAmp),min_H_amp);
                ph_S = RF.protonPhase(j);
                %B1_I = min(round(abs(parameters(j+2))/min_N_amp),1)*sign(parameters(j+2))*B1_C*max(min(abs(parameters(j+2)),N_maxAmp),min_N_amp);
                ph_I = RF.XPhase(j);
                tau  = RF.taus(j);
                
                %%% 1H-Xnuc pulse %%%
                H   = 2*pi*gamma_H*B1_S*(cos(ph_S)*Sx + sin(ph_S)*Sy)+...
                      2*pi*gamma_N*B1_I*(cos(ph_I)*Ix + sin(ph_I)*Iy)+...
                      H_free;
                PH  = kron(eye(2^nspins),H)-kron(H.',eye(2^nspins));
                L   = -1i*PH - diag(RelMtx(:));
                rho = expv(tau,L,rho,1e-7,30);
                
                if (df == 0 && B1_C == 1 && B1_H == 1)
                    rho_ideal = rho;
                end
            end
            rho         = reshape(rho,[size(initState)]);
            WBH         = B1weight_H(bb_H);
            WBC         = B1weight_X(bb_C);
            weight      = WBC*WBH*Wf;
            rem_signal  = trace(reshape(initState,[size(Iz)])'*rho)/traceNorm;
            tred_signal = trace(coil'*rho)/traceNorm; %with the current normalization the equilibrium 1H magnetization is 8
            tred_signal_array(bb_C,bb_H,ff) = tred_signal;
            rem_signal_array(bb_C,bb_H,ff)  = rem_signal;
            sum_weight  = sum_weight + weight;
            waitbar(((ff-1)*lB1r_H*lB1r_X+(bb_C-1)*lB1r_H+bb_H)/(lB0r*lB1r_H*lB1r_X))
            bb_H = bb_H + 1;
        end
        bb_C = bb_C + 1;
    end
end
signal_tred  = sum(sum(sum(tred_signal_array)))/sum_weight;
signal_rem  = sum(sum(sum(rem_signal_array)))/sum_weight;
avgTransfer   = signal_tred
avgPreserved  = signal_rem
close(h);
%% detection
rho           = reshape(rho_ideal,[size(initState)]);
peakTransfer  = trace(coil'*rho)/traceNorm
peakPreserved = trace(Iz'*rho)/traceNorm
nsteps        = 4096;
BW            = 4000;
f             = linspace(-BW/2,BW/2,nsteps);
H_free        = 2*pi*Jcoupling;
PH_free       = kron(eye(2^nspins),H_free)-kron(H_free.',eye(2^nspins));
L_free        = -1i*PH_free - diag(RelMtx(:));
prop_free     = expm(L_free*1/BW);  

for n = 1:nsteps
    fid(n)= trace(coil'*rho)/traceNorm;
    rho   = prop_free*rho(:);
    rho   = reshape(rho,[2^nspins,2^nspins]);
end
%%
fid_apd   = fid.*exp(-0.005*linspace(1,nsteps,nsteps));
spectrum  = ifftshift(ifft(fid_apd));
plot(f,imag(spectrum))
sumIntensity = abs(trapz(spectrum))


















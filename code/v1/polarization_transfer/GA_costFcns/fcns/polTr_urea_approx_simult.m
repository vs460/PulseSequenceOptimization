function[fitness] = polTr_urea_approx_simult(population,pars)

tic
[npop,nvar]    = size(population);
gamma_N        = pars.gX;                 % Xnuc gyromagnetic ratio
gamma_H        = pars.gH;                 % 1H gyromagnetic ratio
J_NH           = -90;                     % coupling in Hz
B1_range_H     = pars.B1_range_H;         % pulse amplitude range for compensation
B1_range_X     = pars.B1_range_X;         % pulse amplitude range for compensation
fRange         = pars.fRange;             % offresonance range for compensation
min_N_amp      = pars.min_X_amp;
min_H_amp      = pars.min_H_amp;
N_maxAmp       = pars.X_maxAmp;
H_maxAmp       = pars.H_maxAmp;
lB1r_H         = length(B1_range_H);
lB1r_X         = length(B1_range_X);
lB0r           = length(fRange);
% B0 inhom distribution
LW = max(fRange)-min(fRange);
lineBroadening = normpdf(fRange,0,LW/2.355);
lineBroadening = lineBroadening/max(lineBroadening);
% B1 inhom distribution
B1weight_H = normpdf(B1_range_H-1,0,1.5/2.355);
B1weight_H = B1weight_H/max(B1weight_H);
B1weight_X = normpdf(B1_range_X-1,0,2/2.355);
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
B = cell((2^nspins)^2,1);
R = zeros((2^nspins)^2,1);
RelMtx = zeros(size(Ix));
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

%% initialization
initState = Sz;
traceNorm = trace(Iz'*Iz);
coil      = Ix-1i*Iy;
%% fitness evaluation
fitness = zeros(npop,1);
for i = 1:npop
    parameters = (population(i,1:nvar));
    taus_tmp   = max(0,parameters(5:5:end));
    B1_S_vec   = min(round(abs(parameters(1:5:end-4))/min_H_amp),1).*sign(parameters(1:5:end-4)).*max(min(abs(parameters(1:5:end-4)),pars.H_maxAmp),min_H_amp)
    B1_I_vec   = min(round(abs(parameters(3:5:end-2))/min_N_amp),1).*sign(parameters(3:5:end-2)).*max(min(abs(parameters(3:5:end-2)),pars.X_maxAmp),min_N_amp)
    B1_S_vec   = smoothCurveToConstraint(gamma_H*B1_S_vec',min(taus_tmp(find(taus_tmp>0))),pars.alpha,pars.beta,1);
    B1_I_vec   = smoothCurveToConstraint(gamma_N*B1_I_vec',min(taus_tmp(find(taus_tmp>0))),pars.alpha,pars.beta,1);

    sum_weight = 0;
    rem_signal_array  = zeros(lB1r,lB1r,lB0r);
    tred_signal_array = zeros(lB1r,lB1r,lB0r);
    for ff = 1:length(fRange)
        df     = fRange(ff);
        H_free = 2*pi*df*offres + 2*pi*Jcoupling;
        bb_C     = 1;
        for B1_C = B1_range_X
            bb_H = 1;
            for B1_H = B1_range_H
                rho  = initState;
                %%%%%%%%%% Pulse sequence %%%%%%%%%%
                for j = 1:5:nvar
                    B1_S = B1_H*B1_S_vec((j-1)/5+1);
                    B1_I = B1_C*B1_I_vec((j-1)/5+1);
                    ph_S = parameters(j+1);
                    ph_I = parameters(j+3);
                    tau  = max(0,parameters(j+4));
                                  
                    H    = 2*pi*B1_I*gamma_N*(cos(ph_I)*Ix + sin(ph_I)*Iy)+...    
                           2*pi*B1_S*gamma_H*(cos(ph_S)*Sx + sin(ph_S)*Sy)+...        
                           H_free;                                                          
                    PH   = expm(-1i*H*tau);
                    PHt  = PH';
                    PrelaxU = ((exp(-RelMtx*tau/2)));
                    rho  = PrelaxU.*rho;
                    rho  = PH*rho*PHt;
                    rho  = PrelaxU.*rho;
                end
                Wf  = lineBroadening(ff);
                WBH = B1weight_H(bb_H);
                WBC = B1weight_X(bb_C);
                weight      = Wf*WBH*WBC;
                tred_signal = trace(coil'*rho)/(trace(Iz'*Iz)); %with the current normalization the equilibrium 1H magnetization is 8             
                tred_signal_array(bb_C,bb_H,ff) = weight*tred_signal;
                sum_weight  = sum_weight + weight;
                bb_H = bb_H + 1;
            end
            bb_C = bb_C + 1;
        end
    end
    signal_tred  = sum(sum(sum((tred_signal_array))))/sum_weight;
    fitness(i)	 = 1e3*gather(-abs(signal_tred));
end
toc
end
















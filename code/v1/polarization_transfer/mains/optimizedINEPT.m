% Genetic algorithm optimized pulse design
% (2019) Vencel Somai -> vs460@cam.ac.uk

clear all
close all
rng default
INEPTcheck = true;
trdir      = 'H2C';
intlChs    = false;       % true = no simultaneous multichannel transmit
%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% !!!!! save the pars structure at the end !!!!! %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gX         = -4.316e6;    % Xnuc gyromagnetic ratio [Hz/T]
gH         = 42.57e6;     % 1H gyromagnetic ratio   [Hz/T]
J          = -90;         % coupling in Hz
T1_I       = 24;          % Xnuc longitudinal
T2_I       = 1.6;         % Xnuc transverse
T1_S       = 2.6;         % 1H longitudinal
T2_S       = 0.06;        % 1H transverse
maxLength  = 2*abs(1/J);  % maximal sequence time length [s]
spinsystem = 'IS';        % 'IS' or 'IS3' and only used for calculating the corresponding INEPT anyway
B1_range_H = 0.2:0.4:2.6; % B1 grid the transfer is evaluated on
B1_range_X = 0.2:0.4:3.4; % B1 grid the transfer is evaluated on
fRange     = -100:25:100; % df grid the transfer is evaluated on
pw         = 1e-4;        % pulse width in the INEPT/HINDER individual 
minPw_I    = 1e-6;        % minimum Xnuc pulse width to achieve frequency selectivity
minPw_S    = 1e-6;        % minimum 1H pulse width to achieve frequency selectivity 
CHch       = 200e-6;      % time for frequency channel change 
% old smoothness parameters, not used anymore
% gwpar      = 10;
% smpar      = 5;
% medfiltpar = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha      = 4e8;         % parameters for smoothness constraint [1st der]
beta       = 4e10;        % parameters for smoothness constraint [2nd der]
%% plot options
plotOpt    = {@gaplotbestf};

%% population options
nvar       = 5*80;
H_maxAmp   = 0.25e-4;
X_maxAmp   = 2e-4;
pop.Type   = 'doubleVector';
pop.Size   = 32;
pop.Create = @gacreationuniform;
% safety check whether the number of different variables are equal
if intlChs
    if mod(nvar,6)~=0
        error('nvar should be multiple of 6')
    end
    pop.Range  = repmat([[0;H_maxAmp],[0;2*pi],[4e-6;maxLength/(nvar/6)],[0;X_maxAmp],[0;2*pi],[4e-6;maxLength/(nvar/6)]],[1,(nvar)/6]);
else
    if mod(nvar,5)~=0
        error('nvar should be multiple of 5')
    end
    pop.Range  = repmat([[0;H_maxAmp],[0;2*pi],[0;X_maxAmp],[0;2*pi],[4e-6;maxLength/(nvar/5)]],[1,(nvar)/5]);
end
population = rand(pop.Size,nvar).*repmat(pop.Range(2,:)-pop.Range(1,:),[pop.Size,1]) + repmat(pop.Range(1,:),[pop.Size,1]);

%% conventional INEPT sequence to validate the simulations
% the last individual in the population is the INEPT sequecne, therefore
% the polarization transfer should be maximal
if INEPTcheck
    testPars = createTestIndividual(nvar,pw,J,gH,gX,spinsystem,intlChs,trdir);
    population(1,:) = testPars;
end

%% selection options
select.Type = @selectionroulette;
ec          = 2;
crossFrac   = 0.8;
%% mutation options
mutFcn      = @mutationuniform;
Scale       = 0.1;
Shrink      = 0.5;
%% crossover options
crssover    = 'crossoverheuristic';
%% number of iterations
NofGens     = 3e5;

%% create options input
if size(population,2) < nvar
    population = [population,zeros(size(population,1),nvar - size(population,2))];
end
pars.gX         = gX;         pars.gH         = gH;
pars.X_maxAmp   = X_maxAmp;   pars.H_maxAmp   = H_maxAmp;
pars.minPw_I    = minPw_I;    pars.minPw_S    = minPw_S;           
pars.X90amp     = 1/4/gX/pw;  pars.H90amp     = 1/4/gH/pw;
pars.min_X_amp  = 4.186e-08;  pars.min_H_amp  = 7.677e-08;
pars.B1_range_H = B1_range_H; pars.fRange     = fRange;
pars.B1_range_X = B1_range_X; pars.CHch       = CHch;
pars.T1_I       = T1_I;       pars.T2_I       = T2_I;
pars.T1_S       = T1_S;       pars.T2_S       = T2_S;
pars.gwpar      = gwpar;      pars.medfiltpar = medfiltpar;
pars.smpar      = smpar;
pars.alpha      = alpha;      pars.beta       = beta;
%%
options = optimoptions('ga','InitialPopulationMatrix',population,'PlotFcn',plotOpt,'PopulationType',pop.Type,'PopulationSize',pop.Size,...
    'CreationFcn',pop.Create,'FitnessScalingFcn',{@fitscalingrank},'InitialPopulationRange',pop.Range,'SelectionFcn',{select.Type},'EliteCount',ec,...
    'CrossoverFraction',crossFrac,'MutationFcn',{mutFcn,Scale},'CrossoverFcn',{crssover},'MaxGenerations',NofGens,...
    'UseVectorized',true,'MaxStallGenerations',NofGens);
%% run optimization and save the solution
if intlChs % no simultaneous pulsing of differetn freq. channels
    parameters = ga(@(x) polTr_urea_approx_intl(x,pars),nvar,[],[],[],[],[],[],[],options);
    idx = 1;
    B1_S_vec = min(round(abs(parameters(1:5:end-4))/pars.min_H_amp),1).*sign(parameters(1:5:end-4)).*max(min(abs(parameters(1:5:end-4)),pars.H_maxAmp),pars.min_H_amp)
    B1_I_vec = min(round(abs(parameters(3:5:end-2))/pars.min_X_amp),1).*sign(parameters(3:5:end-2)).*max(min(abs(parameters(3:5:end-2)),pars.X_maxAmp),pars.min_X_amp)
    %%%%%%%%%%%%% Old smoothness constaint, not used anymore %%%%%%%%%%%%%%
    %     B1_S_vec = medfilt1(smooth(B1_S_vec,pars.smpar),pars.medfiltpar);
    %     B1_I_vec = medfilt1(smooth(B1_I_vec,pars.smpar),pars.medfiltpar);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:6:nvar
        RF.protonAmp(idx)   = B1_S_vec((j-1)/6+1);
        %RF.protonAmp(idx)   = min(round(abs(parameters(j))/pars.min_H_amp),1)*sign(parameters(j))*max(min(abs(parameters(j)),pars.H_maxAmp),pars.min_H_amp);
        RF.protonPhase(idx) = parameters(j+1);
        RF.protonTaus(idx)  = max(0,parameters(j+2));
        RF.XAmp(idx)        = B1_I_vec((j-1)/6+1);
        %RF.XAmp(idx)        = min(round(abs(parameters(j+3))/pars.min_X_amp),1)*sign(parameters(j+3))*max(min(abs(parameters(j+3)),pars.X_maxAmp),pars.min_X_amp);
        RF.XPhase(idx)      = parameters(j+4);
        RF.XTaus(idx)       = min(max(pars.minPw_I,parameters(j+5)),1/abs(J));
        idx = idx+1;
    end
    pars.parameters = parameters;
    pars.RF = RF;
else % simultaneous pulsing on different freq. channels
    parameters = ga(@(x) polTr_urea_approx_simult(x,pars),nvar,[],[],[],[],[],[],[],options);
    %%
    idx = 1;
    taus_tmp = min(max(1e-17,parameters(5:5:end)),1/abs(J));
    T        = cumsum(taus_tmp);
    B1_S_vec = min(round(abs(parameters(1:5:end-4))/pars.min_H_amp),1).*sign(parameters(1:5:end-4)).*max(min(abs(parameters(1:5:end-4)),pars.H_maxAmp),pars.min_H_amp);
    B1_I_vec = min(round(abs(parameters(3:5:end-2))/pars.min_X_amp),1).*sign(parameters(3:5:end-2)).*max(min(abs(parameters(3:5:end-2)),pars.X_maxAmp),pars.min_X_amp);
    B1_S_vec = smoothCurveToConstraint(pars.gH*B1_S_vec',min(taus_tmp(find(taus_tmp>0))),pars.alpha,pars.beta,1);
    B1_I_vec = smoothCurveToConstraint(pars.gX*B1_I_vec',min(taus_tmp(find(taus_tmp>0))),pars.alpha,pars.beta,1);
    
    for j = 1:5:nvar
        %RF.protonAmp(idx)   = min(round(abs(parameters(j))/pars.min_H_amp),1)*sign(parameters(j))*max(min(abs(parameters(j)),pars.H_maxAmp),pars.min_H_amp);
        RF.protonAmp(idx)   = B1_S_vec((j-1)/5+1)/pars.gH;
        RF.protonPhase(idx) = parameters(j+1);
        % RF.XAmp(idx)        = min(round(abs(parameters(j+2))/pars.min_X_amp),1)*sign(parameters(j+2))*max(min(abs(parameters(j+2)),pars.X_maxAmp),pars.min_X_amp);
        RF.XAmp(idx)        = B1_I_vec((j-1)/5+1)/pars.gX;
        RF.XPhase(idx)      = parameters(j+3);
        RF.taus(idx)        = taus_tmp((j-1)/5+1);
        idx = idx+1;
    end
    pars.parameters = parameters;
    pars.RF = RF;
end
%%%%%%%%% !!!save pars structure!!! %%%%%%%%%%%



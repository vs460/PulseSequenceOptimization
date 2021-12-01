% Genetic algorithm optimized pulse design
% (2019) Vencel Somai -> vs460@cam.ac.uk

clear all
close all
rng default
%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% !!!!! save the pars structure at the end !!!!! %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gX         = 10.71e6; %-4.316e6; % Xnuc gyromagnetic ratio [Hz/T]
gH         = 42.57e6;            % 1H gyromagnetic ratio   [Hz/T]
J          = 140;                % coupling in Hz
% relaxation parameters
T1_I       = 7;                  % Xnuc longitudinal
T2_I       = 0.3;                % Xnuc transverse
T1_S       = 1.6;                % 1H longitudinal
T2_S       = 0.1;                % 1H transverse
maxLength  = 1*abs(1/J);         % maximal sequence time length [s]
spinsystem = 'IS';
B1_range_H = 0.5:0.05:1.5;       % B1 grid the transfer is evaluated on
B1_range_X = 0.5:0.05:1.5;       % exp(linspace(-2,2,10));  % B1 grid the transfer is evaluated on
fRange     = -50:10:50;          % df grid the transfer is evaluated on
pw         = 1e-3;               % pulse width in the INEPT/HINDER individual 
minPw_I    = 0e-3;               % minimum Xnuc pulse width to achieve frequency selectivity
minPw_S    = 0;                  % minimum 1H pulse width to achieve frequency selectivity 
alpha      = 8e7;                % parameters for smoothness constraint [1st der]
beta       = 8e7;                % parameters for smoothness constraint [2nd der]
%% plot options
plotOpt    = {@gaplotbestf};

%% population options
nvar       = 150;
H_maxAmp   = 0.25e-4;
X_maxAmp   = 2e-4;
pop.Type   = 'doubleVector';
pop.Size   = 32;
pop.Create = @gacreationuniform;
pop.Range  = repmat([[0;H_maxAmp],[0;2*pi],[4e-6;maxLength/(nvar/3)]],[1,(nvar)/3]);
population = rand(pop.Size,nvar).*repmat(pop.Range(2,:)-pop.Range(1,:),[pop.Size,1]) + repmat(pop.Range(1,:),[pop.Size,1]);
% create selective, off-resonant hard pulse and include it in the initial
% population if required
T = 10; FA = pi/2; offset = 0;
testPars   = createTestJimmuneRF(T,FA,offset,gX,nvar);
population(1,:) = testPars;
%% selection options
select.Type = @selectionstochunif;
ec          = 2;
crossFrac   = 0.8;
%% mutation options
mutFcn      = @mutationuniform;
Scale       = 0.3;
Shrink      = 0.5;
%% crossover options
crssover    = 'crossoverheuristic';
%% number of iterations
NofGens     = 3e4;

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
pars.B1_range_X = B1_range_X; 
pars.T1_I       = T1_I;       pars.T2_I       = T2_I;
pars.T1_S       = T1_S;       pars.T2_S       = T2_S;
pars.alpha      = alpha;      pars.beta       = beta;
%%
options = optimoptions('ga','InitialPopulationMatrix',population,'PlotFcn',plotOpt,'PopulationType',pop.Type,'PopulationSize',pop.Size,...
    'CreationFcn',pop.Create,'FitnessScalingFcn',{@fitscalingrank},'InitialPopulationRange',pop.Range,'SelectionFcn',{select.Type},'EliteCount',ec,...
    'CrossoverFraction',crossFrac,'MutationFcn',{mutFcn,Scale},'CrossoverFcn',{crssover},'MaxGenerations',NofGens,...
    'UseVectorized',true,'MaxStallGenerations',NofGens);
%% run optimization and save the solution
parameters = ga(@(x) JimmuneRF_lac2C_approx(x,pars),nvar,[],[],[],[],[],[],[],options);
%% putting together the result
taus_tmp   = min(max(0,parameters(3:3:end)),1/abs(J));
B1_I_vec   = min(round(abs(parameters(1:3:end-2))/pars.min_X_amp),1).*sign(parameters(1:3:end-2)).*max(min(abs(parameters(1:3:end-2)),pars.X_maxAmp),pars.min_X_amp)
B1_I_vec   = smoothCurveToConstraint(gX*B1_I_vec',min(taus_tmp(find(taus_tmp>0))),pars.alpha,pars.beta,1);
idx        = 1;
% and storing it in the pars structure
for j = 1:3:nvar-2
    RF.Amp(idx)   = B1_I_vec((j-1)/3+1)/gX;
    RF.Phase(idx) = parameters(j+1);
    RF.Taus(idx)  = min(max(0,parameters(j+2)),1/abs(J));
    idx = idx+1;
end
pars.parameters = parameters;
pars.RF = RF;
% !!!!! save pars !!!!!
%% simulate result if required
testJimmuneRF_ureaN15(pars)

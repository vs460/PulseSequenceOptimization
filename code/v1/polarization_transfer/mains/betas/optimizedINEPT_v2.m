% Genetic algorithm optimized pulse design
%
% This version also adds an optional frequency offset to the pulses
%
% (2020) Vencel Somai -> vs460@cam.ac.uk

clear all
close all
rng default
INEPTcheck = true;
globalSearch = false;
%% parameters
gamma_C = 10.71e6;              % 13C gyromagnetic ratio
gamma_H = 42.57e6;              % 1H gyromagnetic ratio
J = 4.1;                        % coupling in Hz
maxLength = 1/J;                % maximal sequence time length  in s
interlievedChannels = true;     % true = Varian design(no simultaneous pulsing)
spinsystem = 'I3S';
%% plot options
plotOpt = {@gaplotbestf};

%% population options
nvar       = 240
% maximal rotation/subpulse to constrain the maximal pulse amplitude
H_maxAmp   = 0.25e-4;
C_maxAmp   = 1e-4;
pop.Type   = 'doubleVector';
pop.Size   = 32;
pop.Create = @gacreationuniform;
% safety check whether the number of different variables are equal
if interlievedChannels
    if mod(nvar,8)~=0
        error('nvar should be multiple of 6')
    end
    pop.Range  = repmat([[0;H_maxAmp],[0;2*pi],[4e-6;maxLength/((nvar)/6)],[-2000;2000],[0;C_maxAmp],[0;2*pi],[4e-6;maxLength/((nvar)/8)],[-2000;2000]],[1,(nvar)/8]);
else
    if mod(nvar,5)~=0
        error('nvar should be multiple of 5')
    end
    pop.Range  = repmat([[0;H_maxAmp],[0;2*pi],[0;C_maxAmp],[0;2*pi],[4e-6;maxLength/((nvar)/5)]],[1,(nvar)/5]);
end
population = rand(pop.Size,nvar).*repmat(pop.Range(2,:)-pop.Range(1,:),[pop.Size,1]) + repmat(pop.Range(1,:),[pop.Size,1]);
%population = [];
%% conventional INEPT sequence to validate the simulations
% the last individual in the population is the INEPT sequecne, therefore
% the polarization transfer should be maximal
if INEPTcheck
    if interlievedChannels
        pulseTime = 1e-3;
        if strcmp(spinsystem,'IS')
            tauMix = 1/4/J-pulseTime;
            tauPrep =  1/4/J-pulseTime;
        elseif strcmp(spinsystem,'I3S')
            tauPrep = 1/4/J;
            tauMix = 0.5*acos(sqrt(2/3))/pi/J;
        else
            error('unknown spin system');
        end
        INEPTpopulation = zeros(1,nvar);
        INEPTpopulation(1:8) = [0,0,0,0,1/4/gamma_C/pulseTime,0,pulseTime,0];                % first 90x pulse on proton
        INEPTpopulation(9:16) = [0,0,0,0,0,0,tauMix,0];  % tau/2 free evolution
        INEPTpopulation(17:24) = [2/4/gamma_H/pulseTime,0,pulseTime,0,2/4/gamma_C/pulseTime,0,pulseTime,0];   % 180x pulse on both spins
        INEPTpopulation(25:32) =  [0,0,0,0,0,0,tauMix,0];   % tau/2 free evolution
        INEPTpopulation(33:40) = [1/4/gamma_H/pulseTime,0,pulseTime,0,1/4/gamma_C/pulseTime,pi/2,pulseTime,0];  % second 90y pulse on both spins
        INEPTpopulation(41:48) = [0,0,0,0,0,0,tauPrep,0];                      % tau'/2 free evolution
        INEPTpopulation(49:56) = [2/4/gamma_H/pulseTime,0,pulseTime,0,2/4/gamma_C/pulseTime,0,pulseTime,0];   % second 180x on both spins
        INEPTpopulation(57:64) = [0,0,0,0,0,0,tauPrep,0];                  % tau'/2 free evolution
        INEPTpopulation(65:end) = 0;
%         HINDERpopulation = zeros(1,nvar);
%         HINDERpopulation(1:6) = [1/4/gamma_H/pulseTime,pi/2,pulseTime,1/4/gamma_C/pulseTime,pi/2,pulseTime];% first 90x pulse on proton
%         HINDERpopulation(7:12) = [0,0,tauMix/2,0,0,tauMix/2];  % tau/2 free evolution         7*pi/4 fro half preserved, -0.058*pi fro 0.13 preserved
%         HINDERpopulation(13:18) = [2/4/gamma_H/pulseTime,3*pi/2,pulseTime,2/4/gamma_C/pulseTime,7*pi/4,pulseTime];% 180x pulse on both spins
%         HINDERpopulation(19:24) =  [0,0,tauMix/2,0,0,tauMix/2];% tau/2 free evolution
%         HINDERpopulation(25:30) = [1/4/gamma_H/pulseTime,pi/2,pulseTime,1/4/gamma_C/pulseTime,pi/2,pulseTime];% second 90y pulse on both spins
%         HINDERpopulation(31:36) = [0,0,tauPrep,0,0,tauPrep]; % tau'/2 free evolution
%         
%         HINDERpopulation(37:42) = [0,0,0,1/4/gamma_C/pulseTime,pi/2,pulseTime];% second 180x on both spins
%         HINDERpopulation(43:48) = [0,0,0,2/4/gamma_C/pulseTime,0,pulseTime];% second 180x on both spins
%         HINDERpopulation(49:54) = [0,0,0,1/4/gamma_C/pulseTime,pi/2,pulseTime];% second 180x on both spins
%         
%         HINDERpopulation(55:60) = [1/4/gamma_H/pulseTime,pi/2,pulseTime,0,0,0];% second 180x on both spins
%         HINDERpopulation(61:66) = [2/4/gamma_H/pulseTime,0,pulseTime,0,0,0];% second 180x on both spins
%         HINDERpopulation(67:72) = [1/4/gamma_H/pulseTime,pi/2,pulseTime,0,0,0];% second 180x on both spins
%        
%         HINDERpopulation(73:78) = [0,0,tauPrep,0,0,tauPrep]; % tau'/2 free evolution
%         HINDERpopulation(79:end) = 0;
         population(end,:) = INEPTpopulation;
        %population = repmat(HINDEpopulation,[pop.Size,1]);
    else
        if strcmp(spinsystem,'IS')
            tauMix = 1/4/J;
        elseif strcmp(spinsystem,'I3S')
            tauMix = 0.5*acos(sqrt(2/3))/pi/J;
        else
            error('unknown spin system');
        end
        pulseTime = 1e-3;
        INEPTpopulation = zeros(1,nvar);
        INEPTpopulation(1:5) = [1/4/gamma_C/pulseTime,0,pulseTime,0,0,0];   % first 90x pulse on proton
        INEPTpopulation(6:10) = [0,0,0,0,tauMix];  % tau/2 free evolution
        INEPTpopulation(11:15) = [1/4/gamma_H/pulseTime,0,1/4/gamma_C/pulseTime,0,2*pulseTime];   % 180x pulse on both spins
        INEPTpopulation(16:20) =  [0,0,0,0,tauMix];   % tau/2 free evolution
        INEPTpopulation(21:25) = [1/4/gamma_H/pulseTime,pi/2,1/4/gamma_C/pulseTime,0,pulseTime];  % second 90y pulse on both spins
        INEPTpopulation(26:30) = [0,0,0,0,1/4/J];                      % tau'/2 free evolution
        INEPTpopulation(31:35) = [1/4/gamma_H/pulseTime,0,1/4/gamma_C/pulseTime,0,2*pulseTime];   % second 180x on both spins
        INEPTpopulation(36:40) = [0,0,0,0,1/4/J];                  % tau'/2 free evolution
        INEPTpopulation(41:end) = 0;
        population(end,:) = INEPTpopulation;
        population = repmat(INEPTpopulation,[pop.Size,1]);
    end
end
%population_new(1:size(population,1),:) = population;
%population = population_new;
%% selection options
select.Type = @selectionroulette;
ec          = 2;
crossFrac   = 0.8;
%% mutation options
mutFcn     = @mutationgaussian;
Scale      = 0.50;
Shrink     = 0.5;
%% crossover options
crssover    = 'crossoverintermediate';
%% number of iterations
NofGens     = 3e4;

%% create options input
%'InitialPopulationMatrix',repmat(pulse_phase,[pop.Size,1]),
%load('population3'); %   population = [population;rand(pop.Size,nvar).*repmat(pop.Range(2,:),[pop.Size,1])]; pop.Size = pop.Size*2;
if size(population,2) < nvar
    population = [population,zeros(size(population,1),nvar - size(population,2))];
end
pars.C_maxAmp  = C_maxAmp;
pars.H_maxAmp  = H_maxAmp;
pars.C90amp    = 1/4/gamma_C/pulseTime;
pars.H90amp    = 1/4/gamma_H/pulseTime;
pars.min_C_amp = 0.2*pars.C90amp;
pars.min_H_amp = 0.2*pars.H90amp;
%%
options = optimoptions('ga','InitialPopulationMatrix',population,'PlotFcn',plotOpt,'PopulationType',pop.Type,'PopulationSize',pop.Size,...
    'CreationFcn',pop.Create,'InitialPopulationRange',pop.Range,'SelectionFcn',{select.Type},'EliteCount',ec,...
    'CrossoverFraction',crossFrac,'MutationFcn',{mutFcn,Scale,Shrink},'CrossoverFcn',{crssover},'MaxGenerations',NofGens,...
    'UseVectorized',true,'MaxStallGenerations',NofGens);
%,'HybridFcn','fminsearch'
%% run optimization and save the solution
if interlievedChannels
    parameters = ga(@(x) polTransfed_lactate_approx_intl_C2H_v2(x,pars),nvar,[],[],[],[],[],[],[],options);
    %parameters = ga(@(x) polarizationTransfered_lactate_relax_interleaved(x,pars),nvar,[],[],[],[],[],[],[],options);
    RF.carbonAmp = parameters(1:6:end-5);
    RF.carbonPhase = parameters(2:6:end-4)
    RF.carbonTaus = parameters(3:6:end-3)
    RF.protonAmp = parameters(4:6:end-2)
    RF.protonPhase = parameters(5:6:end-1)
    RF.protonTaus = parameters(6:6:end)
    save(strcat('C:\Users\Somai01\Documents\pulseSequenceDesign\INEPT_optimization\INEPT_rf_interleaved'));
else
    %parameters = ga(@polarizationTransfered_universal_relax_diffBasis,nvar,[],[],[],[],[],[],[],options);
    %parameters = ga(@polarizationTransfered_universal_relax,nvar,[],[],[],[],[],[],[],options);
    %parameters = ga(@polarizationTransfered_lactate_relax,nvar,[],[],[],[],[],[],[],options);
    parameters = ga(@(x) polTransfed_lactate_approximate(x,pars),nvar,[],[],[],[],[],[],[],options);
    RF.protonAmp = parameters(1:5:end-4);
    RF.protonPhase = parameters(2:5:end-3)
    RF.carbonAmp = parameters(3:5:end-2)
    RF.carbonPhase = parameters(4:5:end-1)
    RF.taus = parameters(5:5:end)
    save(fullfile(pwd,'INEPT_rf.mat'),RF);
end



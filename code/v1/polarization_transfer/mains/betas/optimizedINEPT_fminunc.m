% Genetic algorithm optimized pulse design
% (2019) Vencel Somai -> vs460@cam.ac.uk

clear all
close all
rng default
INEPTcheck   = true;
globalSearch = false;
%% parameters
gX   = 10.71e6;          % Xnuc gyromagnetic ratio [Hz/T]
gH   = 42.57e6;          % 1H gyromagnetic ratio   [Hz/T]
J    = 4.1;              % coupling in Hz
% relaxation parameters
T1_I = 50;               % Xnuc longitudinal
T2_I = 0.3;              % Xnuc transverse
T1_S = 1.73;             % 1H longitudinal
T2_S = 0.1;              % 1H transverse
maxLength  = abs(1/J);   % maximal sequence time length [s]
intlChs    = true;       % true = no simultaneous multichannel transmit
spinsystem = 'I3S';
B1_range   = 0.1:0.3:2;  % B1 grid the transfer is evaluated on
fRange     = -200:40:200;% df grid the transfer is evaluated on
pw         = 0.5e-3;     % pulse width in the INEPT/HINDER individual 
%% plot options
plotOpt    = {@gaplotbestf};

%% population options
nvar       = 60;
H_maxAmp   = 0.25e-4;
X_maxAmp   = 1e-4;
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
    testPars = createTestIndividual(nvar,pw,J,gH,gX,spinsystem,intlChs);
    population(end,:) = testPars;
end

%% selection options
select.Type = @selectionroulette;
ec          = 2;
crossFrac   = 0.8;
%% mutation options
mutFcn      = @mutationgaussian;
Scale       = 0.5;
Shrink      = 0.5;
%% crossover options
crssover    = 'crossoverheuristic';
%% number of iterations
NofGens     = 3e4;

%% create options input
if size(population,2) < nvar
    population = [population,zeros(size(population,1),nvar - size(population,2))];
end
pars.gX        = gX;       pars.gH        = gH;
pars.X_maxAmp  = X_maxAmp; pars.H_maxAmp  = H_maxAmp;
pars.X90amp    = 1/4/gX/pw;pars.H90amp    = 1/4/gH/pw;
pars.min_X_amp = 4.186e-08;pars.min_H_amp = 7.677e-08;
pars.B1_range  = B1_range; pars.fRange    = fRange;
pars.T1_I      = T1_I;     pars.T2_I      = T2_I;
pars.T1_S      = T1_S;     pars.T2_S      = T2_S;

%% run the optimization
options = optimoptions('fminunc','Display','iter','MaxFunEvals',1e8,'MaxIter',5e1,'TolFun',1e-10,'TolX',1e-10,'PlotFcn','optimplotfval','SpecifyObjectiveGradient',false);%'SpecifyObjectiveGradient',true
pulse   = fminunc(@(pulse) polTransfed_lactate_approximate_fminunc(pulse,pars),INEPTpopulation,options);

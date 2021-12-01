% Script to interpolate the designed pulse so that every slot is integer
% multiple of a given timestep set by the user. Finer the timestep the
% less modification is made in the pulse sequence.

clear all
close all
%% load the optimization result
[file,path] = uigetfile({'*.*'});
load(fullfile(path,file));
%%
RF = pars.RF;
nvar        = 80;                 % number of slots in the pulse, 1/5 or 1/6 of the nvar in optimizedINEPT script
fileWrite   = false;               % set true for write the interpolated pulses into a Bruker compatibel textfile
NaCalib = false;
%% parameters
gamma_X = 10.71e6;                % 13C gyromagnetic ratio in Hz
gamma_H = 42.57e6;                % 1H gyromagnetic ratio in Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pw_Na23    = 104e-6;               %%%%%%%%%%%%%%%%% Set %%%%%%%%%%%%%%%%%%
Na23offset = -3830;             %%%%%%%%%%%%%%%% these %%%%%%%%%%%%%%%%%
WaterOffset= -825;              %%%%%%%%%%%%%%%%% !!! %%%%%%%%%%%%%%%%%%
% calibrated parameters for Xnuc
P_X        = 45;                  % pulse power [dB]
T_X        = 284e-6;               % corresponding 90 degree pulse length [s] lactate
% calibrated parameters for proton
P_H        = 45;                  % pulse power [dB]
T_H        = 3*515e-6;              % corresponding 90 degree pulse length [s] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eqPw_C13   = 68.8/62.8*Pw_Na23

if NaCalib
    T_X = eqPw_C13;
    C13offset_Na = (79.282e6+Na23offset)/(79.282e6-3555)*(75.369e6-2014)-75.369e6
    C13offset_H  = (299.709e6+WaterOffset)/(299.709e6-1084)*(75.369e6-2285)-75.369e6
end
fieldStrength_X = 1/4/T_X;
fieldStrength_H = 1/4/T_H;

%% proton channel
TS = 4e-6;                    % sampling timestep of the interpolated pulses
H_amp_interpolated   = [];
H_phase_interpolated = []; 

% looping through the slots
for i = 1:nvar
    timeSlot_amp = RF.protonAmp(i);
    if timeSlot_amp < 0                                                 % if the amplitude is negative change the phase with pi
        timeSlot_phase = mod(RF.protonPhase(i)+pi,2*pi);
        timeSlot_amp   = abs(timeSlot_amp);
    else
        timeSlot_phase = RF.protonPhase(i);
    end
    N_of_points        = round(RF.taus(i)/TS);                                 % number of points needed for the upsampling
    ampScaleFactor     = (RF.taus(i)/TS)/round(RF.taus(i)/TS);              % scaling factor to correct the amplitude of the point so that the integral remains same
    timeSlot_amp       = repmat(timeSlot_amp*ampScaleFactor,[1,N_of_points]); % correct the amplitudes with the scaling factor and split the slot to N_of_points timestep
    timeSlot_phase     = repmat(timeSlot_phase,[1,N_of_points]);            % same for the phase
    H_amp_interpolated = [H_amp_interpolated,timeSlot_amp];             % storing the interpolated points to the final pulse variable
    H_phase_interpolated = [H_phase_interpolated,timeSlot_phase];
end
% writing to file according to Bruker format (amplitude = percentage, phase = [0,360] in degrees)
if fileWrite
    fileID = fopen('protonRF','w');
    for n = 1:length(H_amp_interpolated)
        A = H_amp_interpolated(n)*100/max(abs(H_amp_interpolated));
        P = H_phase_interpolated(n)*180/pi;
        fprintf(fileID,'%5.2f, %5.2f\n',A,P);
    end
    fclose(fileID);
    % write shape to file with vnmrj file format
    vNMRjwrite(H_amp_interpolated.*exp(1i*H_phase_interpolated),1,'vs_C2lacH2C_Hshape_INEPT.RF');
end
%% carbon channel
X_amp_interpolated   = [];
X_phase_interpolated = [];
for i = 1:nvar
    timeSlot_amp = RF.XAmp(i);
    if timeSlot_amp < 0
        timeSlot_phase = mod(RF.XPhase(i)+pi,2*pi);
        timeSlot_amp   = abs(timeSlot_amp);
    else
        timeSlot_phase = RF.XPhase(i);
    end
    N_of_points        = round(RF.taus(i)/TS);
    ampScaleFactor     = (RF.taus(i)/TS)/round(RF.taus(i)/TS);
    timeSlot_amp       = repmat(timeSlot_amp*ampScaleFactor,[1,N_of_points]);
    timeSlot_phase     = repmat(timeSlot_phase,[1,N_of_points]);
    X_amp_interpolated = [X_amp_interpolated,timeSlot_amp];
    X_phase_interpolated = [X_phase_interpolated,timeSlot_phase];
end

if fileWrite
    fileID = fopen('nitorgenRF','w');
    for n = 1:length(X_amp_interpolated)
        A = X_amp_interpolated(n)*100/max(abs(X_amp_interpolated));
        P = X_phase_interpolated(n)*180/pi;
        fprintf(fileID,'%5.2f, %5.2f\n',A,P);
    end
    fclose(fileID);
    % write shape to file with vnmrj file format
    vNMRjwrite(X_amp_interpolated.*exp(1i*X_phase_interpolated),1,'vs_C2lacH2C_Cshape_INEPT.RF');
end
%% display parameters
seqtime = TS*length(H_amp_interpolated)     % total time of the sequence in S
max(abs(H_amp_interpolated))*gamma_H;       % proton channel maximal amplitude in Hz
max(abs(X_amp_interpolated))*gamma_X;       % carbon channel maximal amplitude in Hz
protonAmp_dB    = P_H + 20*log10(max(abs(RF.protonAmp))*gamma_H/fieldStrength_H)
nitroAmp_dB     = P_X + 20*log10(abs(max(abs(RF.XAmp))*gamma_X/fieldStrength_X))

%% plot pulse
% alpha = 4e6;beta = 4e6;
% RF.protonAmp = 1e-4*smoothCurveToConstraint(RF.protonAmp'*1e4,min(RF.taus(find(RF.taus>0))),alpha,beta,1)';
% RF.XAmp = 1e-4*smoothCurveToConstraint(RF.XAmp'*1e4,min(RF.taus(find(RF.taus>0))),alpha,beta,1)';
 
figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%% RF H chanel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1)
xlabel('Time [ms]')
title('Optimized RF - ^1H CH')
yyaxis left
plot(1e3*cumsum(RF.taus),1e4*abs(RF.protonAmp),'k','LineWidth',1.5);pbaspect([1,1,1])
%plot(linspace(0,1e3*sum(RF.Taus),length(X_amp_interpolated)),X_amp_interpolated)
ylabel('Pulse amplitude [G]')
ylim([0,1.2*max(abs(RF.protonAmp)*1e4)])
set(gca,'FontWeight','bold','YColor','k')
yyaxis right
AX=plot(1e3*cumsum(RF.taus),(RF.protonPhase).*((-1).^(RF.protonAmp<0)),'--','Color',[.5 .5 .5],'LineWidth',1.5);pbaspect([1,1,1])
%plot(linspace(0,1e3*sum(RF.Taus),length(X_amp_interpolated)),X_phase_interpolated,'--')
ylabel('Pulse phase [rad]')
set(gca,'FontWeight','bold','YColor',[.5 .5 .5]) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% RF X chanel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2)
xlabel('Time [ms]')
title('Optimized RF - ^{13}C CH')
yyaxis left
plot(1e3*cumsum(RF.taus),1e4*abs(RF.XAmp),'k','LineWidth',1.5);pbaspect([1,1,1])
%plot(linspace(0,1e3*sum(RF.Taus),length(X_amp_interpolated)),X_amp_interpolated)
ylabel('Pulse amplitude [G]')
ylim([0,1.2*max(abs(RF.XAmp)*1e4)])
set(gca,'FontWeight','bold','YColor','k')
yyaxis right
AX=plot(1e3*cumsum(RF.taus),(RF.XPhase).*((-1).^(RF.XAmp<0)),'--','Color',[.5 .5 .5],'LineWidth',1.5);pbaspect([1,1,1])
%plot(linspace(0,1e3*sum(RF.Taus),length(X_amp_interpolated)),X_phase_interpolated,'--')
ylabel('Pulse phase [rad]')
set(gca,'FontWeight','bold','YColor',[.5 .5 .5])
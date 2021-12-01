% Script to interpolate the designed pulse so that every slot is integer
% multiple of a given timestep set by the user. Finer the timestep the
% less degradation is made in the pulse sequence by the interpolation.

clear all
close all
%% load the optimization result (the pars.mat file)
[file,path] = uigetfile({'*.*'});
load(fullfile(path,file));
%%
NaCalib    = false;   % whether the body sodium is used for centre freq. and Tx gain calibration
RF         = pars.RF;
nvar       = 50;      % number of slots/pulse points in the pulse, 1/3 of the nvar number
shape_name = 'vs_C2lacJimmuneRF9p76ms.RF';
fileWrite  = false;   % set true for write the interpolated pulses into VnmrJ textfile
%% parameters
gamma_X    = 10.71e6; % 13C gyromagnetic ratio in Hz
% calibrated parameters for Xnuc
Pw_Na23    = 80e-6;   % 90deg pulse width for Na23
Na23offset = -3867.6; % Centre freq. of body sodium
WaterOffset= -896.8;  % Centre freq. of body water
eqPw_C13   = 68.8/62.8*Pw_Na23 % equivalent 13C hard pulse length
% if Na23 calibration is not used measure these
P_X        = 45;      % pulse power [dB]
T_X        = 335e-6;  % corresponding 90 degree pulse length [s]
if NaCalib
    T_X = eqPw_C13;
    C13offset_Na = (79.282e6+Na23offset)/(79.282e6-3555)*(75.369e6-2014)-75.369e6
    C13offset_H = (299.709e6+WaterOffset)/(299.709e6-1084)*(75.369e6-2285)-75.369e6
end
fieldStrength_X = 1/4/T_X;

TS = 4e-6;            % sampling timestep of the interpolated pulses
%% carbon channel
X_amp_interpolated     = [];
X_phase_interpolated   = [];
for i = 1:nvar          % interpolate each slot/pulse point to TS resolution
    timeSlot_amp       = RF.Amp(i);
    if timeSlot_amp < 0  % negate the phase if the amplitude is negative
        timeSlot_phase = mod(RF.Phase(i)+pi,2*pi);
        timeSlot_amp   = abs(timeSlot_amp);
    else
        timeSlot_phase = RF.Phase(i);
    end
    N_of_points        = round(RF.Taus(i)/TS);
    % correct for interpolation error
    ampScaleFactor     = (RF.Taus(i)/TS)/round(RF.Taus(i)/TS);
    timeSlot_amp       = repmat(timeSlot_amp*ampScaleFactor,[1,N_of_points]);
    timeSlot_phase     = repmat(timeSlot_phase,[1,N_of_points]);
    % concatenate the result
    X_amp_interpolated = [X_amp_interpolated,timeSlot_amp];
    X_phase_interpolated = [X_phase_interpolated,timeSlot_phase];
end

if fileWrite
    % write shape to file with vnmrj file format
    vNMRjwrite(X_amp_interpolated.*exp(1i*X_phase_interpolated),1,shape_name);
end
%% display parameters
seqtime = TS*length(X_amp_interpolated) % total time of the sequence in S
Amp_dB  = P_X + 20*log10(abs(max(abs(RF.Amp))*gamma_X/fieldStrength_X)) % Tx gain [dB]
%% plot final shape
figure
xlabel('Time [ms]')
title('Optimized RF')
%%%
yyaxis left
hold on
%bar(1e3*cumsum(RF.Taus+1e-10),1e4*abs(RF.Amp))
plot(linspace(TS,1e3*sum(RF.Taus),length(X_amp_interpolated)),1e4*X_amp_interpolated,'b')
ylabel('Pulse amplitude [G]')
ylim([0,1.2*max(X_amp_interpolated*1e4)])
set(gca,'FontWeight','bold')
%%%
yyaxis right
%plot(1e3*cumsum(RF.Taus),(RF.Phase).*((-1).^(RF.Amp<0)))
plot(linspace(0,1e3*sum(RF.Taus),length(X_amp_interpolated)),X_phase_interpolated,'--')
ylabel('Pulse phase [rad]')
%%%
set(gca,'FontWeight','bold')
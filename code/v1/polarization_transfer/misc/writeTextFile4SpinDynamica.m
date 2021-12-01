% Writes the pulses into text file for SpinDynamica validation

clear all; close all
%% load pulse sequence 
[file,path] = uigetfile({'*.*'});
load(fullfile(path,file))
fileWrite   = false;       
spinsystem  = 'C2lac';

%% parameters
parameters  = pars.parameters;
RF          = pars.RF;
nvar        = length(parameters)/5;
gamma_X     = pars.gX;                 
gamma_H     = pars.gH;                 
reH         = real(RF.protonAmp.*exp(1i*RF.protonPhase))*gamma_H;
imH         = imag(RF.protonAmp.*exp(1i*RF.protonPhase))*gamma_H;
reX         = real(RF.XAmp.*exp(1i*RF.XPhase))*gamma_X;
imX         = imag(RF.XAmp.*exp(1i*RF.XPhase))*gamma_X;

%% write file
path   = uigetdir();
%%
fileID = fopen(fullfile(path,spinsystem),'w');
fprintf(fileID,'%s = ','taus');
fprintf(fileID,'%6.8f,',RF.taus);
fprintf(fileID,'\n');
fprintf(fileID,'%s = ','reH');
fprintf(fileID,'%6.2f,',reH);
fprintf(fileID,'\n');
fprintf(fileID,'%s = ','imH');
fprintf(fileID,'%6.2f,',imH);
fprintf(fileID,'\n');
fprintf(fileID,'%s = ','reX');
fprintf(fileID,'%6.2f,',reX);
fprintf(fileID,'\n');
fprintf(fileID,'%s = ','imX');
fprintf(fileID,'%6.2f,',imX);
fclose(fileID);



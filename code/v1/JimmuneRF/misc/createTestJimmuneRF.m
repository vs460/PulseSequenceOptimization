function [testPars] = createTestJimmuneRF(T,FA,offset,gamma,nvar)
% helper function to create selective, off-resonant hard pulse if needed
testPars            = zeros(1,nvar);
T                   = T/1e3;             % conversion from [s] to [ms]
dt                  = T/(nvar/3); 
A                   = FA/(2*pi)/gamma/T; % conversio from flip-angle to [T]
phase               = 2*pi*offset*linspace(0,T,nvar/3);
testPars(1:3:end-2) = A;                 % Amplitude
testPars(2:3:end-1) = phase;             % phase
testPars(3:3:end)   = dt;                % time steps
end
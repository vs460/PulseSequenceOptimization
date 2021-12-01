function [testPars] = createTestIndividual(nvar,pw,J,gH,gX,spinsystem,intlChs,trdir)
J = abs(J);
if intlChs
    if strcmp(spinsystem,'IS')
        tauMix  = 1/4/J-pw;
        tauPrep =  1/4/J-pw;
    elseif strcmp(spinsystem,'I3S')
        tauMix  = 1/4/J;
        tauPrep = 0.5*acos(sqrt(2/3))/pi/J;
    else
        error('unknown spin system');
    end
    INEPTpop         = zeros(1,nvar);
    if strcmp(trdir,'C2H')
    INEPTpop(1:6)    = [0,0,0,1/4/gX/pw,0,pw];                % first 90x pulse on proton
    elseif strcmp(trdir,'H2C')
            INEPTpop(1:6)    = [1/4/gH/pw,0,pw,0,0,0];                % first 90x pulse on proton
    else
        error('Undefined transfer direction')
    end
    INEPTpop(7:12)   = [0,0,0,0,0,tauMix];  % tau/2 free evolution
    INEPTpop(13:18)  = [2/4/gH/pw,0,pw,2/4/gX/pw,0,pw];   % 180x pulse on both spins
    INEPTpop(19:24)  = [0,0,0,0,0,tauMix];   % tau/2 free evolution
    if strcmp(trdir,'C2H')
        INEPTpop(25:30)  = [1/4/gH/pw,0,pw,1/4/gX/pw,pi/2,pw];  % second 90y pulse on both spins
    elseif strcmp(trdir,'H2C')
        INEPTpop(25:30)  = [1/4/gH/pw,pi/2,pw,1/4/gX/pw,0,pw];  % second 90y pulse on both spins
    else
        error('Undefined transfer direction')
    end
    INEPTpop(31:36)  = [0,0,0,0,0,tauPrep];                      % tau'/2 free evolution
    INEPTpop(37:42)  = [2/4/gH/pw,0,pw,2/4/gX/pw,0,pw];   % second 180x on both spins
    INEPTpop(43:48)  = [0,0,0,0,0,tauPrep];                  % tau'/2 free evolution
    INEPTpop(49:end) = 0;
    %         HINDERpop         = zeros(1,nvar);
    %         HINDERpop(1:6)    = [1/4/gH/pw,pi/2,pw,1/4/gX/pw,pi/2,pw];% first 90x pulse on proton
    %         HINDERpop(7:12)   = [0,0,tauMix/2,0,0,tauMix/2];  % tau/2 free evolution         7*pi/4 fro half preserved, -0.058*pi fro 0.13 preserved
    %         HINDERpop(13:18)  = [2/4/gH/pw,3*pi/2,pw,2/4/gX/pw,7*pi/4,pw];% 180x pulse on both spins
    %         HINDERpop(19:24)  = [0,0,tauMix/2,0,0,tauMix/2];% tau/2 free evolution
    %         HINDERpop(25:30)  = [1/4/gH/pw,pi/2,pw,1/4/gX/pw,pi/2,pw];% second 90y pulse on both spins
    %         HINDERpop(31:36)  = [0,0,tauPrep,0,0,tauPrep]; % tau'/2 free evolution
    %
    %         HINDERpop(37:42)  = [0,0,0,1/4/gX/pw,pi/2,pw];% second 180x on both spins
    %         HINDERpop(43:48)  = [0,0,0,2/4/gX/pw,0,pw];% second 180x on both spins
    %         HINDERpop(49:54)  = [0,0,0,1/4/gX/pw,pi/2,pw];% second 180x on both spins
    %
    %         HINDERpop(55:60)  = [1/4/gH/pw,pi/2,pw,0,0,0];% second 180x on both spins
    %         HINDERpop(61:66)  = [2/4/gH/pw,0,pw,0,0,0];% second 180x on both spins
    %         HINDERpop(67:72)  = [1/4/gH/pw,pi/2,pw,0,0,0];% second 180x on both spins
    
    %         HINDERpop(73:78)  = [0,0,tauPrep,0,0,tauPrep]; % tau'/2 free evolution
    %         HINDERpop(79:end) = 0;
    testPars = INEPTpop;
    %population = repmat(HINDEpopulation,[pop.Size,1]);
else
    if strcmp(spinsystem,'IS')
        tauPrep = 1/4/J; %391e-6;%1/4/J;
        tauMix  = 1/4/J;
        delta   = 18.050/180*pi;
    elseif strcmp(spinsystem,'I3S')
        tauPrep = 0.5*acos(sqrt(2/3))/pi/J;
        tauMix  =  1/4/J;
    else
        error('unknown spin system');
    end
    INEPTpop         = zeros(1,nvar);
    INEPTpop(1:5)    = [1/4/gH/pw,0,0,0,pw];   % first 90x pulse on proton
    INEPTpop(6:10)   = [0,0,0,0,tauMix-3/2*pw];  % tau/2 free evolution
    INEPTpop(11:15)  = [1/4/gH/pw,0,1/4/gX/pw,0,2*pw];   % 180x pulse on both spins
    INEPTpop(16:20)  = [0,0,0,0,tauMix-3/2*pw];   % tau/2 free evolution
    INEPTpop(21:25)  = [1/4/gH/pw,pi/2,1/4/gX/pw,0,pw];  % second 90y pulse on both spins
    INEPTpop(26:30)  = [0,0,0,0,tauPrep-3/2*pw];                      % tau'/2 free evolution
    INEPTpop(31:35)  = [1/4/gH/pw,0,1/4/gX/pw,0,2*pw];   % second 180x on both spins
    INEPTpop(36:40)  = [0,0,0,0,tauPrep-pw];                  % tau'/2 free evolution
    INEPTpop(41:end) = 0;
    testPars = INEPTpop;
    
    %         HINDERpop         = zeros(1,nvar);
    %         HINDERpop(1:5)    = [1/4/gH/pw,pi/2,  1/4/gN/pw,pi/2,        pw];% first 90x pulse on proton
    %         HINDERpop(6:10)   = [0,        0,     0,        0,           tauMix];  % tau/2 free evolution         7*pi/4 fro half preserved, -0.058*pi fro 0.13 preserved
    %         HINDERpop(11:15)  = [2/4/gH/pw,3*pi/2,2/4/gN/pw,3*pi/2+delta,pw];% 180x pulse on both spins
    %         HINDERpop(16:20)  = [0,        0,     0,        0,           tauMix];% tau/2 free evolution
    %         HINDERpop(21:25)  = [1/4/gH/pw,pi/2,  1/4/gN/pw,pi/2,        pw];% second 90y pulse on both spins
    %         HINDERpop(26:30)  = [0,        0,     0,        0,           tauPrep]; % tau'/2 free evolution
    %         HINDERpop(31:35)  = [1/4/gH/pw,pi/2,  1/4/gN/pw,pi/2,        pw];% second 180x on both spins
    %         HINDERpop(36:40)  = [2/4/gH/pw,0,     2/4/gN/pw,0,           pw];% second 180x on both spins
    %         HINDERpop(41:45)  = [1/4/gH/pw,pi/2,  1/4/gN/pw,pi/2,        pw];% second 180x on both spins
    %         HINDERpop(46:50)  = [0,        0,     0,        0,           tauPrep]; % tau'/2 free evolution
    %         HINDERpop(79:end) = 0;
    %         testPars = HINDERpop;
end


end


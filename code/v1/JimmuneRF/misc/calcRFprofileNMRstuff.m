function [cost,mz_array,freqOffset] = calcRFprofileNMRstuff(B1_I_vec,parameters,Taus,pars,eval)
% Runs Blcoh simulation if frequency selectivity is required, e.g. in an
% hyperpolarized in vivo lactate detection experiment 


%% simulate pulse
% putting together the pulse and scaling it to units of G
pulse_tmp = 1e4*B1_I_vec/pars.gX.*exp(1i*parameters(2:3:end-1))';
npoints   = length(pulse_tmp);
% if the function used for evaluating an optimization result set fine grid
if eval == 1
    freqOffset      = linspace(-1000,12000,13000);
    pars.B1_range_X = linspace(0,max(pars.B1_range_X),200);
else % if used in the optimization use coarser grid for speed
    freqOffset      = linspace(-1000,12000,1300);
end
    
goal      = ones(1,length(freqOffset)); weight = goal-1;
weight(find(abs(freqOffset-10245)<200)) =1;

ii = 1;
% evaluation on B1 grid
for B1 = pars.B1_range_X
    %% initial conditions
    mx = zeros(1,length(freqOffset));my = zeros(1,length(freqOffset));mz = ones(1,length(freqOffset));
    pulse = B1*pulse_tmp;
    [mx_tmp,my_tmp,mz_tmp] = bloch1(pulse,0,Taus,pars.T1_I,pars.T2_I,freqOffset,0,[0 0]);
    while sum(isnan(mx_tmp)+isnan(my_tmp)+isnan(mz_tmp))>0
        [mx_tmp,my_tmp,mz_tmp] = bloch1(pulse,0,Taus,pars.T1_I,pars.T2_I,freqOffset,0,[0 0],mx,my,mz);
    end
    mx = mx_tmp; my = my_tmp; mz = mz_tmp;
    mx_array(ii,:) = mx;
    my_array(ii,:) = my;
    mz_array(ii,:) = mz;
    ii = ii+1;
end
%% calculate cost function as the squared error between goal and the actual freq. prof.
cost = sum(sum(weight.*(mz_array-repmat(goal,[length(pars.B1_range_X),1])).^2))/sum(weight);

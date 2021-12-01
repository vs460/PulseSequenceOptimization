%% simulate pulse
addpath(your bloch simulator...
pulse_tmp = B1_I_vec/gamma_C.*exp(1i*parameters(2:3:end-1))';
npoints   = length(pulse_tmp);
Taus      = min(max(0,parameters(3:3:end)),1/abs(J_CH));

fRange    = linspace(-1000,12000,130);
for B1 = B1_range_X
    %% initial conditions
    mx = zeros(1,length(fRange));my = zeros(1,length(fRange));mz = ones(1,length(fRange));
    pulse = 1i*B1*pulse_tmp;
    [mx_tmp,my_tmp,mz_tmp] = bloch1(pulse,0,Taus,T1_I,T2_I,fRange,0,[0 0],mx,my,mz);
    while sum(isnan(mx_tmp)+isnan(my_tmp)+isnan(mz_tmp))>0
        [mx_tmp,my_tmp,mz_tmp] = bloch1(pulse,0,Taus,T1_I,T2_I,fRange,0,[0 0],mx,my,mz);
    end
    mx = mx_tmp; my = my_tmp; mz = mz_tmp;
    mx_array(i,:) = mx;
    my_array(i,:) = my;
    mz_array(i,:) = mz;
    i = i+1;
end

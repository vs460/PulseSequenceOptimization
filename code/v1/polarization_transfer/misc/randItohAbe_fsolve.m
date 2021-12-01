function [y] = randItohAbe_fsolve(x,d,tau_min,tau_max,pars)
% Derivative free solver, only beta version

tau = sqrt(tau_min*tau_max);
%% define the operators
%%%%%%%%%%%% This is the cost function %%%%%%%%%%%%%%%
V   = @(x) polTransfed_urea_approx_N2H(x,pars);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F   = @(beta) beta + (V(x-tau*beta*d)-V(x))/(tau*beta);

%% try to solve for the stepsize
try
    beta = fsolve(F,10);
catch 
    y = x;
    return;
end
y  = x - tau*beta*d;
Vx = V(x);
Vy = V(y);
while (abs(Vy-Vx)/norm(y-x).^2 > 1/tau_min || abs(Vy-Vx)/norm(y-x).^2 < 1/tau_max)
    if abs(Vy-Vx)/norm(y-x).^2 > 1/tau_min
        y = y/sigma;
    else
        y = y*sigma;
    end
end
end

        
         

    


  




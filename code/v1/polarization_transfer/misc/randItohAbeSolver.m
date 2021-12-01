function [y] = randItohAbeSolver(x,d,tau_min,tau_max,pars)
% Derivative free solver, only beta
sigma = 0.9;
tau = sqrt(tau_min*tau_max);
eps = 1e-2;
V = @(x) polTransfed_urea_approx_N2H(x,pars);
Vx = V(x);
Vxepsd = V(x+eps*d);
if Vxepsd>=Vx
    d = -d;
    Vxepsd = V(x+eps*d);
    if Vxepsd>=Vx
        y = x;
        return;
    end
end

beta = -(Vxepsd-Vx)/(eps*tau);
x0 = x; x1 = x+eps*d; x2 = x+beta*d;
Vx0 = V(x0); Vx1 = V(x1);Vx2 = V(x2);
while (Vx2-Vx1)./norm(x2-x1) <= (Vx1-Vx0)./norm(x1-x0)
    beta = beta/sigma;
    x2 = x+beta*d;
    Vx2 = V(x2);
end
y = x1;
y = x1 - 0.5*((x1-x0).^2*(Vx1-Vx2)-(x1-x2).^2*(Vx1-Vx0))./...
             ((x1-x0)*(Vx1-Vx2)-(x1-x2)*(Vx1-Vx0));
y(find(isnan(y))) = 0;
while V(y)>= Vx
    if V(y) >= V(x1)
        x0 = x1; x1 = y; x2 = x2;
    else
        x0 = x0; x2 = x1; x1 = y;
    end
    Vx0 = V(x0); Vx1 = V(x1);Vx2 = V(x2);
    y = x1 - 0.5*((x1-x0).^2*(Vx1-Vx2)-(x1-x2).^2*(Vx1-Vx0))/...
             ((x1-x0)*(Vx1-Vx2)-(x1-x2)*(Vx1-Vx0));
    y(find(isnan(y))) = 0;
end
Vy = V(y);
while (abs(Vy-Vx)/norm(y-x)^2 > 1/tau_min || abs(Vy-Vx)/norm(y-x)^2 < 1/tau_max)
    if abs(Vy-Vx)/norm(y-x)^2 > 1/tau_min
        y = y/sigma;
    else
        y = y*sigma;
    end
end
V(y)
end

        
         

    


  


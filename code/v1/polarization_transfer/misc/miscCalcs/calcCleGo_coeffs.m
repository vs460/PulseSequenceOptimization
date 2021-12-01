function [result] = CG_coeff(j1,m1,j2,m2,j,m)
% Calculates Clebsch-Gordan coefficients

if j2 ~=1
    error("single spin rank must be 1")
end

switch m2
    case -1
        if j == j1 + 1
            if m1 == m + 1
                coeff = sqrt((j1-m)*(j1-m+1)/(2*j1+1)/(2*j1+2));
            else
                error("k numbers do not match")
            end
        elseif j == j1 
            if m1 == m + 1
                coeff = sqrt((j1-m)*(j1+m+1)/(2*j1)/(j1+1));
            else
                error("k numbers do not match")
            end
        elseif j == j1 - 1 
            if m1 == m + 1
                coeff = sqrt((j1+m+1)*(j1+m)/(2*j1)/(2*j1+1));
            else
                error("k numbers do not match")
            end
        else
            error("no such case for q and k")
        end
    case 0
        if j == j1 + 1
            if m1 == m
                coeff = sqrt((j1-m+1)*(j1+m+1)/(2*j1+1)/(j1+1));
            else
                error("k numbers do not match")
            end
        elseif j == j1 
            if m1 == m
                coeff = m/sqrt(j1*(j1+1));
            else
                error("k numbers do not match")
            end
        elseif j == j1 - 1 
            if m1 == m
                coeff = -sqrt((j1-m)*(j1+m)/(j1)/(2*j1+1));
            else
                error("k numbers do not match")
            end
        else
            error("no such case for q and k")
        end
    case 1
        if j == j1 + 1
            if m1 == m - 1 
                coeff = sqrt((j1+m)*(j1+m+1)/(2*j1+1)/(2*j1+2));
            else
                error("k numbers do not match")
            end
        elseif j == j1 
            if m1 == m - 1
                coeff = -sqrt((j1+m)*(j1-m+1)/(2*j1)/(j1+1));
            else
                error("k numbers do not match")
            end
        elseif j == j1 - 1 
            if m1 == m - 1
                coeff = sqrt((j1-m)*(j1-m+1)/(2*j1)/(2*j1+1));
            else
                error("k numbers do not match")
            end
        else
            error("no such case for q and k")
        end
    otherwise
        error("wrong single spin q number")
end
result = coeff;
end

        


function [relax_Ham] = calc_Relax_op_lactate(R,B)

% Calculates the basis coefficients of the relaxtion superoperator R
% in the given superoperator basis specified by B. Then the effective
% Hamiltonian superoperator can be back-converted to the Hilber-space to
% reduce the computational burden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Change the path for saving the results %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(B)       % number of basis elements
coeff = zeros(N,1);
relax_Ham = zeros(sqrt(N),sqrt(N));
h = waitbar(0,'Calculating operator coefficients...');
for i = 1:N
    op_tmp = cell2mat(B(i));
    sup_op_tmp = zeros(N,N);
    parfor r = 1:N
        for s = 1:N
            sup_op_tmp(r,s) = trace(cell2mat(B(r))'*(op_tmp*cell2mat(B(s))-cell2mat(B(s))*op_tmp));
        end
    end
    if (trace(sup_op_tmp'*sup_op_tmp) == 0)
        coeff(i) = 0;
    else
        coeff(i) = trace(R'*sup_op_tmp)/trace(sup_op_tmp'*sup_op_tmp);
    end
    relax_Ham = relax_Ham + coeff(i)*cell2mat(B(i));
    save(fullfile('C:\Users\Somai01\Documents\pulseSequenceDesign\INEPT_optimization\lactate\sup_operator_basis',strcat('basis_element_',num2str(i))),'sup_op_tmp');
    waitbar(i/(N*N),h,['Calculating operator coefficients: ' sprintf('%i%% along',round(i/(N*N)*100))]);
end
close(h)
save('C:\Users\Somai01\Documents\pulseSequenceDesign\INEPT_optimization\lactate\sup_operator_basis\basis_coeffs.mat','coeff');
save('C:\Users\Somai01\Documents\pulseSequenceDesign\INEPT_optimization\lactate\sup_operator_basis\basis_set.mat','B');
save('C:\Users\Somai01\Documents\pulseSequenceDesign\INEPT_optimization\lactate\sup_operator_basis\relax_Ham.mat','relax_Ham');



end


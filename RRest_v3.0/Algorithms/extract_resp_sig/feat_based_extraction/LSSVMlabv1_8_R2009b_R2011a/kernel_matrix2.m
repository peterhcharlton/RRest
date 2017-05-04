function omega = kernel_matrix2(omega,kernel_type, kernel_pars,d)
% INTERNAL USE ONLY
% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

if strcmp(kernel_type,'RBF_kernel')
    omega = exp(-omega./(2*kernel_pars)); 
  
elseif strcmp(kernel_type,'RBF4_kernel')
    omega = 0.5*(3-omega./kernel_pars).*exp(-omega./(2*kernel_pars));
    
% elseif strcmp(kernel_type,'sinc_kernel')
%     omega = sinc(omega./kernel_pars); 
    
% elseif strcmp(kernel_type,'wav_kernel')  
%     omega = cos(kernel_pars(3)*omega{2}./kernel_pars(2)).*exp(-omega{1}./kernel_pars(1)); 

elseif strcmp(kernel_type,'lin_kernel')
    return

elseif strcmp(kernel_type,'poly_kernel')
    omega = (omega + kernel_pars(1)).^kernel_pars(2);
end
    
  
function [irf,vardecomp] = SVAR(Y,p,constant,irf_transform,shock_index,H)

% Author: Alvin Lo Hei Chun 

% Given the identified B matrix, compute the IRF and Variance Decomposition

% Using Bootstrap to compute the Confidence bound of IRF

% If H is not specified, H is computed using Cholesky Decomposition

%==========================================================================

[T,N] = size(Y);

var_result = VAR2(Y,p,constant);

%--------------------------------------------------
% Impulse Response Function
%--------------------------------------------------
T_irf = 20;
shock = H(:,shock_index);
irf = SVAR_bootstrap(Y,p,shock,irf_transform);

%--------------------------------------------------
% Variance Decomposition
%--------------------------------------------------
psi_H = zeros(N,N,T_irf);
for ii = 1:T_irf
    psi_H(:,:,ii) = var_result.vma(:,:,ii)*H;
end
psi = squeexe(psi_H(:,shock_index,:));
variance = H(:,shock_index)*H(:,shock_index)';
MSE = var_result.variance;
vardecomp = zeros(N,T_irf-1);

for jj = 1:T_irf-1
    vardecomp(:,jj) = diag(variance).*(diag(MSE)^(-1));
    MSE = MSE + psi_H(:,:,jj+1)*psi_H(:,:,jj+1)';
    variance = variance + psi(:,jj+1)*psi(:,jj+1)';
end
end




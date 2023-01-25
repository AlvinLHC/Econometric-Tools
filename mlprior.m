function [beta_prior,Sigma_prior] = mlprior(Y,p,gamma,w) 

% Author: Alvin Lo Hei Chun 
% Compute the Litterman Prior 
if nargin == 2
    gamma = 0.2;
    w = 0.5;
end

N = size(Y,2);

% Compute tau as suggested by Hamilton (1994)
tau = zeros(N,1);
for n = 1:N 
    single_var  = arfit(Y(:,n),p,1); % With Constant
    tau(n)      = sqrt(single_var.variance);
end

Sigma_prior = zeros(N,N,p);

for i_p = 1:p
    for i_n = 1:N
        for j_n = 1:N 
            if i_n == j_n 
                Sigma_prior(i_n,j_n,i_p) = gamma^2/i_p^2;
            else
                Sigma_prior(i_n,j_n,i_p) = w*gamma*tau(i_n)/(tau(j_n)*i_p);
            end
        end
    end
end

beta_prior = zeros(N,N,p);
beta_prior(:,:,1) = eye(N);

% Construct the stacked form
beta_prior = reshape(beta_prior,N,N*p);
beta_prior = reshape(beta_prior',N*N*p,1);

Sigma_prior_diag = reshape(Sigma_prior,N,N*p);
Sigma_prior_diag = reshape(Sigma_prior_diag',N*N*p,1);
Sigma_prior = diag(Sigma_prior_diag);
end
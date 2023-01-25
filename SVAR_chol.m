function [irf_result,H] = SVAR_chol(Y,p,shock_index,constant,irf_transform)

% Author: Alvin Lo Hei Chun 

% Using Cholesky Decomposition to do the structural VAR
% Plot IRF given a shock to the variable in column (index)

% Confidence Interval of IRF:
%   1) Monte-Carlo
%   2) Bootstrap

% INPUT: 
%==========================================================================

[T,N] = size(Y);

if nargin == 3
    constant = 0;
end

var_result = VAR2(Y,p,constant);
sigma_U = var_result.variance;
H = chol(sigma_U,'lower');
shock = H(:,shock_index);

% Impulse Response Function
%--------------------------------------------------
T_irf = 20;
irf_pre = zeros(N,T_irf);
for t = 1:T_irf 
    irf_pre(:,t) = var_result.vma(:,:,t)*shock;
end

%==========================================================================
% Monte-Carlo Impulse Response
%==========================================================================
X = lag_form(Y,p);
Phi = var_result.beta_companion;
if constant == 1
    X = [ones(T-p,1), X];
    Phi = [var_result.constant,Phi(1:N,:)]';
else
    Phi = Phi(1:N,:)';
end
phi = Phi(:);
Q = X'*X/(T-p);
sigma_phi = kron(sigma_U,inv(Q))/(T-p); % variance of vec(Phi)

Dn = dupmat(N);
Dn_plus = inv(Dn'*Dn)*Dn';
sigma_omega = 2*Dn_plus*kron(sigma_U,sigma_U)*Dn_plus'/T; % variance of vech(omega)

% Ensure the positive definitness of sigma_Phi and sigma_omega by replacing
% the negative eigenvalues by zero
[~,pd_Phi] = cholcov(sigma_phi);
while pd_Phi ~=0
    [P_phi,V_phi] = eig(sigma_phi);
    V_phi = max(V_phi,0);
    sigma_phi = P_phi*V_phi/P_phi;
    [~,pd_Phi] = cholcov(sigma_phi);
end

[~,pd_omega] = cholcov(sigma_omega);
while pd_omega ~=0
    [P_omega,V_omega] = eig(sigma_omega);
    V_omega = max(V_omega,0);
    sigma_omega = P_omega*V_omega/P_omega;
    [~,pd_omega] = cholcov(sigma_omega);
end

% Drawing 
monte_carlo_draw = 10000;
monte_carlo_phi = mvnrnd(phi,sigma_phi,monte_carlo_draw)';
monte_carlo_omega = mvnrnd(vech(sigma_U),sigma_omega,monte_carlo_draw)';
irf_store = zeros(monte_carlo_draw,T_irf,N);
% Given drawings, compute impulse response function 
%---------------------------------------------------
% (Draw of omega might not be positive definite, use the technique above to
% transform it into positive definite) 

for i_mc = 1:monte_carlo_draw
    phi_draw = monte_carlo_phi(:,i_mc);
    omega_draw = monte_carlo_omega(:,i_mc);
    omega_draw = reshape(Dn*omega_draw,N,N);
    irf_draw = irf_chol(phi_draw,omega_draw,p,constant,shock_index);
    irf_store(i_mc,:,:) = irf_draw;
end

% Transform the IRF
first_diff = find(irf_transform == 1);
irf_store(:,:,first_diff) = cumsum(irf_store(:,:,first_diff),2);
second_diff = find(irf_transform == 2);
irf_store(:,:,second_diff) = cumsum(cumsum(irf_store(:,:,second_diff),2),2);

irf_result = zeros(T_irf,3,N);
quant = [0.1,0.5,0.9];
for i_n = 1:N
    irf_result(:,:,i_n) = quantile(irf_store(:,:,i_n),quant)';
end

% % Transform the IRF
% first_diff = find(irf_transform == 1);
% irf_result(:,:,first_diff) = cumsum(irf_result(:,:,first_diff));
% second_diff = find(irf_transform == 2);
% irf_result(:,:,second_diff) = cumsum(cumsum(irf_result(:,:,second_diff)));
end

function series = irf_chol(beta,sigma_U,p,constant,shock_index)

% Given VAR coefficients and sigma_U, compute the impulse response function
% using Cholesky Decomposition

N = size(sigma_U,2);
[~,pd_U] = cholcov(sigma_U);
while pd_U ~=0
    [P_u,V_u] = eig(sigma_U);
    V_u = max(V_u,0);
    sigma_U = P_u*V_u/P_u;
    [~,pd_U] = cholcov(sigma_U);
end
H = chol(sigma_U,'lower');
shock = H(:,shock_index);
Phi = reshape(beta,N*p+constant,N)';
if constant == 1
    phi0 = Phi(:,1);
    Phi = Phi(:,2:end);
end
Phi_companion = [Phi;eye(N*(p-1)),zeros(N*(p-1),N)];

T_irf = 20;
series = zeros(N,T_irf);

for t = 1:T_irf 
    temp = Phi_companion^(t-1);
    series(:,t) =temp(1:N,1:N)*shock;
end
series = series';
end


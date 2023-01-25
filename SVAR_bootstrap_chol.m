function irf_bootstrap_result = SVAR_bootstrap_chol(Y,p,shock_index,irf_transform)

B = 10000; % Number of times to do bootstrap
% con_bound = 0.025*B ;

[T,N] = size(Y);
T_irf = 20;

var_result = VAR2(Y,p,1);

U = var_result.residual;
U_index = randi([1,T-p],T+40,B);
tau_index = randi([p,T],B,1);
irf_bootstrap = zeros(B,T_irf,N);

for b = 1:B
    e = U(U_index(:,b),:);
    Ysim = varsim(var_result,tau_index(b),e');
    var_result_bootstrap = VAR2(Ysim,p,1);
    beta = var_result_bootstrap.beta_stack;
    sigma_U = var_result_bootstrap.variance;
    irf_bootstrap(b,:,:) = irf_chol(beta,sigma_U,p,1,shock_index);
    if mod(b,1000) == 0
        disp(['Bootstrap ', int2str(b)]);
    end
end

% Transform IRF
%------------------------------------
first_diff = find(irf_transform == 1);
irf_bootstrap(:,:,first_diff) = cumsum(irf_bootstrap(:,:,first_diff),2);
second_diff = find(irf_transform == 2);
irf_bootstrap(:,:,second_diff) = cumsum(cumsum(irf_bootstrap(:,:,second_diff),2),2);

irf_bootstrap_result = zeros(T_irf,3,N);
quant = [0.1,0.5,0.9];
for n = 1:N
    irf_bootstrap_result(:,:,n) = quantile(irf_bootstrap(:,:,n),quant)';
end

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
        
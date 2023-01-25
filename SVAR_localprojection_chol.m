function irf_bootstrap_result = SVAR_localprojection_chol(Y,p,shock_index,irf_transform)

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
    sigma_U = var_result_bootstrap.variance;
    H = chol(sigma_U,'lower');
    shock = H(:,shock_index);
    irf_sim = localprojection(Ysim,p,shock);
    irf_bootstrap(b,:,:) = irf_sim;
    if mod(b,1000) == 0
        disp(['Bootstrap Local Projection ',int2str(b)]);
    end
end

% Transform IRF
%------------------------------------
first_diff = find(irf_transform == 1);
irf_bootstrap(:,:,first_diff) = cumsum(irf_bootstrap(:,:,first_diff),2);
second_diff = find(irf_transform == 2);
irf_bootstrap(:,:,second_diff) = cumsum(cumsum(irf_bootstrap(:,:,second_diff),2));

irf_bootstrap_result = zeros(T_irf,3,N);
quant = [0.1,0.5,0.9];
for n = 1:N
    irf_bootstrap_result(:,:,n) = quantile(irf_bootstrap(:,:,n),quant)';
end

end

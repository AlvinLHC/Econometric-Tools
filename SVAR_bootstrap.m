function irf_bootstrap_result = SVAR_bootstrap(Y,p,shock,irf_transform)

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
    T_irf = 20;
    irf_pre = zeros(N,T_irf);
    %----------------------------------------------------------------------
    % Impulse Response Function
    %----------------------------------------------------------------------
    for t = 1:T_irf 
        irf_pre(:,t) = var_result_bootstrap.vma(:,:,t)*shock;
    end
    %----------------------------------------------------------------------
    % Variance Decomposition
    %----------------------------------------------------------------------
    % variance = shock*shock';
    
    irf_bootstrap(b,:,:) = irf_pre';
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

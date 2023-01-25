function var_result = VAR2(Y,p,constant) 

% Author: Alvin Lo Hei Chun

% Compute the VAR Result using the tensor product method 
% Y = X beta + e 
% X = korn(I,x)

[T,N] = size(Y);

if nargin == 2
    constant = 0;
end

if constant == 0
    X_reg = lag_form(Y,p);
else 
    X_reg = [ones(T-p,1),lag_form(Y,p)];
end

% Computing the coefficients
%--------------------------------------------------------------------------
Y_reg = Y(p+1:T,:); % (T-p) x N
Y_reg = Y_reg(:); % N(T-p)
X_reg = kron(eye(N,N),X_reg); % (N(T-p))x(N(1+Np))
beta_stack = (X_reg'*X_reg)\X_reg'*Y_reg;  % N(1+Np)x1

% Computing error and variance 
%--------------------------------------------------------------------------
predicted_error_stack = Y_reg - X_reg*beta_stack;
predicted_error = reshape(predicted_error_stack,T-p,N);
sigma_y = predicted_error'*predicted_error/(T-p);

% Reshape the stacked coefficient vector
%--------------------------------------------------------------------------
if constant ~=0
    beta_companion_pre = reshape(beta_stack,N*p+1,N)'; % N x (1+Np)
    beta_companion = [beta_companion_pre(:,2:end);eye(N*(p-1)),zeros(N*(p-1),N)];
    beta_constant = beta_companion_pre(:,1);
else 
    beta_companion = [reshape(beta_stack,N*p,N)';eye((p-1)*N),zeros((p-1)*N,N)];
end

% Compute the VMA matrices
%--------------------------------------------------------------------------
T_s = 20;
vma = zeros(N,N,T_s);

for t = 1:20
    temp = beta_companion^(t-1);
    vma(:,:,t) = temp(1:N,1:N);
end

% Information Criteria 
%--------------------------------------------------------------------------
var_result.IC.AIC = log(det(sigma_y)) + 2*N^2*p/(T-p);
var_result.IC.HQIC = log(det(sigma_y)) + 2*N^2*p*log(log(T-p))./(T-p);
var_result.IC.BIC = log(det(sigma_y)) + N^2*p*log(T-p)/(T-p);

% Check the stability 
%--------------------------------------------------------------------------
stability = max(abs(eig(beta_companion)));

% Store Results
%--------------------------------------------------------------------------
% Store Basic
var_result.basic.T = T;
var_result.basic.p = p;
var_result.basic.Y = Y;
var_result.basic.Y_reg = Y_reg;
var_result.basic.X_reg = X_reg;

% Store result
var_result.beta_stack = beta_stack;
var_result.beta_companion = beta_companion;
if constant == 1
    var_result.constant = beta_constant;
end
var_result.stability = stability;
var_result.residual = predicted_error;
var_result.variance = sigma_y;
var_result.vma = vma;
%--------------------------------------------------------------------------
end



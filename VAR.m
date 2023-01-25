function var_result = VAR(X,p,constant)

% Author: Alvin Lo Hei Chun
% Subject to changes, as this program does not give a constant term

% Coefficient
% Fitted value
% Prediction
% Error
% Covariance Matrix 
% R square 
% Adjusted R square 
% IRF 
% Confidence interval
% Variance of coefficients 
% p value 

if nargin == 2
    constant = 0;
end

[T,N] = size(X);

Y_reg = X((p+1):end,:);

if constant == 0
    X_reg = lag_form(X,p); % No Constant 
else 
    X_reg = [ones(T-p,1),lag_form(X,p)]; 
end

A = zeros(N,N,p);
A_mat = inv(X_reg'*X_reg)'*(X_reg'*Y_reg);

if constant == 0
    A_rest = A_mat;
    A0 = zeros(N,1);
else 
    A0 = A_mat(:,1);
    A_rest = A_mat(2:end,:);
end

for jj = 1:p
    A(:,:,jj) = A_rest(((jj-1)*N+1):jj*N,:)' ;
end

% Three forms of coefficient matricies 
%   1. Compact form: N x N x p + N array
%   2. Extensive form: (1 + Nxp)xN array
%   3. Vectorize form: vectorize (2)

result.p = p;
result.X = X;

result.coeff_compact = A;

result.coeff_extensive = A_mat;

result.coeff_vectorize = A_mat(:);

result.coeff_companion = [A_mat';eye(N*(p-1)),zeros(N*(p-1),N)];

[~,eigval] = eig(result.coeff_companion);
result.stability = max(abs(diag(eigval))) < 1;

result.error = Y_reg - X_reg*A_mat ; 

result.covariance = result.error'*result.error/(T-p);

result.predict = X_reg*A_mat;
result.fitted = reshape([X(p+1:end,:);result.predict],T-p,2,N);

var_result = result;

end


    





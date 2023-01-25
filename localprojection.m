function irf = localprojection(Y,p,d,Ts)

% Use local projection to plot the impulse response function according to
% Oscar Jorda AER (2005)
% Author: Alvin Lo Hei Chun
%
% INPUTS: 
%       1) X: The data series 
%       2) p: The lag length 
%       3) d: The shock
%       4) Ts: IRF plot length, = 20 is not specified 
if nargin < 4
    Ts = 20;
end

[T,N] = size(Y);
irf = zeros(N,Ts);

for s = 1:Ts
    
% Two ways to do the same VAR estimation, uncomment the way you prefer
%--------------------------------------------------------------------------
%     Y = X(p+s:T,:); % Size = (T-p-s+1) x N 
%     % Construct Xs 
%     Xs = zeros(T-s-p+1,N*p);
%     for i_p = 1:p 
%         Xs(:,N*(i_p-1)+1:N*i_p) = X(p+1-i_p:T-s+1-i_p,:);
%     end
%     Xs = [ones(T-p-s+1,1),Xs];
%     A = (Xs'*Xs)\(Xs'*Y);
%     % Size = (Np) x (N)
%     irf(:,s) = (A(2:N+1,:)')*d; % Size = Nx1
%--------------------------------------------------------------------------
    X_reg = [ones(T-p-s+1,1),lag_form(Y(1:T-s+1,:),p)];
    Y_reg = Y(p+s:T,:);
    Y_reg = Y_reg(:);
    X_reg = kron(eye(N),X_reg);
    beta_stack = (X_reg'*X_reg)\X_reg'*Y_reg;
    beta_mat = reshape(beta_stack,N*p+1,N)';
    irf(:,s) = beta_mat(:,2:N+1)*d;
%--------------------------------------------------------------------------
end
irf = irf';

% Monte Carlo Drawing 


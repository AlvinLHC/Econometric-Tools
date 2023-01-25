function series = varsim(var_result,tau,e,Tsim)

% Author: Alvin Lo Hei Chun
% INPUT: 
%   1) var_result:      generated from using VAR2
%   2) T:               Simulation Period, equal to T if not specified
%   3) tau:             Starting period. equal to p if not specified
%   4) e:               Shocks 

p = var_result.basic.p;
Y = var_result.basic.Y;
[T,N] = size(Y);

if nargin < 4
    Tsim = T + 40;
elseif nargin == 1
    Tsim = T + 40;
    tau = p;
end

if nargin < 3
    e = chol(var_result.variance,'lower')*randn(N,Tsim);
end

y0 = Y(tau-p+1:tau,:);
y0 = reshape(y0',N*p,1);
F = var_result.beta_companion;
c = var_result.constant;

ysim = zeros(N,Tsim);
for t = 1:Tsim 
    temp = F*y0;
    ysim(:,t) = temp(1:N) + e(:,t);
    y0 = [ysim(:,t);y0(1:N*(p-1))];
end
series = ysim(:,41:end)';
end

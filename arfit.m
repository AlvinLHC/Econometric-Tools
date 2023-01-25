function result = arfit(y,p,constant)

% Using AR(p) model to estimate the process y 
% y is assumed to be 1D
if nargin == 2
    constant = 0;
end

if p == 0 
    result.coef = 0;
    result.coef_companion = 0;
    result.error = y;
    result.variance = var(y);
else
    T = size(y,1);
    Y = y(p+1:T);
    X = zeros(T-p,p);
    X(:,1) = y(p:T-1);

    for i = 1:p-1
        X(:,i+1) = y(p-i:T-1-i);
    end 
    
    if constant ~= 0
        X = [ones(T-p,1),X];
    end
    
    result.coef = inv(X'*X)*X'*Y;
    if constant ~= 0
        result.constant = result.coef(1);
        result.coef_companion = [result.coef(2:end)';eye(p-1),zeros(p-1,1)];
    else
        result.coef_companion = [result.coef';eye(p-1),zeros(p-1,1)];
    end
    result.error = Y - X*result.coef;
    result.variance = result.error'*result.error/(T-p);

end

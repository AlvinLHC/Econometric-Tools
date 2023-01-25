function X_lag = lag_form(X,p)

% Express a time series in information form with p lags

[T,N] = size(X);

X_lag = zeros(T-p, N*p);

for j = 1:p
    
    X_lag(:,(N*(j-1)+1):N*j) = X((p+1-j):(T-j),:) ;
    
end

end

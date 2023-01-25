function tran_data = data_transformer(X,index1,index2)

% index 1: determine the transformation method
%   1. Level: nothing is applied
%   2. First Difference
%   3. Second Difference
%   4. log
%   5. log first difference
%   6. log second difference


% index2 = 1: hpfilter 

if index1 == 1
    X_data = X;
elseif index1 == 2
    X_data = X(2:end,:) - X(1:end-1,:);
elseif index1 == 3
    X_data = X(3:end,:) - 2* X(2:end-1,:) + X(1:end-2,:);
elseif index1 == 4
    X_data = log(X);
elseif index1 == 5
    X_data = 100*(log(X(2:end,:)) - log(X(1:end-1,:)));
else X_data = log(X(3:end,:)) - 2*log(X(2:end-1,:)) + log(X(1:end-2,:));
end
    
if index2 == 1
    X_data = X_data - hpfilter(X_data, 1600); % default is quarterly data
end
    
T = size(X_data,1);
tran_data = X_data;
end
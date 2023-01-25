% Data transformer 2

function tran_data = data_transformer2(X,tran_matrix)

% Dimension of the tran_matrix should be 2 x N

[T,N] = size(X);
N1 = size(tran_matrix,2);

if N ~= N1
    error('The number of column of X and tran_matrix are not equal')
end

X_data = zeros(T,N);
index = 0;

for j = 1:N
    
    X_pre = data_transformer(X(:,j),tran_matrix(1,j),tran_matrix(2,j));
    
    if T - length(X_pre)== 2
        X_data(2:end-1,j) = X_pre;
        index = max(index,2);
    elseif T - length(X_pre)== 1
        X_data(2:end,j) = X_pre;   
        index = max(index,1);
    else
        X_data(:,j) = X_pre;
    end
end
if index == 2
    X_data = X_data(2:end-1,:);
elseif index == 1
    X_data = X_data(2:end,:);
else
end

tran_data = X_data;

end

    
        

    

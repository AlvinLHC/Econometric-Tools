function data = defaultdata()

% Load default data

quarterly = xlsread('Data','Quarterly','B2:C189');
monthly = xlsread('Data','Monthly','B2:B565');

% Transform monthly to quarterly 
monthly = squeeze(mean(reshape(monthly,3,size(monthly,1)/3,size(monthly,2))));
gdp = 400*diff(log(quarterly(:,1)));
inflation = 400*diff(log(quarterly(:,2)));
%inflation = diff(log(monthly(:,2)));
FFR = monthly';
data = [gdp,inflation,FFR(2:end)];
end



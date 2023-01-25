function [index,info] = varIC(Y)

pmax = 8;
p_grid = (1:pmax)';
AIC = zeros(pmax,1);
HQIC = zeros(pmax,1);
BIC = zeros(pmax,1);

for i_p = 1:pmax
    p = p_grid(i_p);
    var_result = VAR2(Y,p,1);
    AIC(i_p) = var_result.IC.AIC;
    HQIC(i_p) = var_result.IC.HQIC;
    BIC(i_p) = var_result.IC.BIC;
end

[~, a_index] = min(AIC);
[~, h_index] = min(HQIC);
[~, b_index] = min(BIC);
index = [a_index,h_index,b_index]';
info = [AIC,HQIC,BIC];
    
end

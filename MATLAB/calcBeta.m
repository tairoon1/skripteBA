function [beta] = calcBeta(ah,beta_0,ah_min,beta_min,beta_07)

A = [0 0 1
     ah_min^2 ah_min 1
     0.7^2 0.7 1];
b = [beta_0;beta_min;beta_07];
parameters = A\b

beta = parameters(1)*ah.^2 + parameters(2)*ah + parameters(3);


%% plot fitted function
% f = @(x) parameters(1)*x.^2 + parameters(2)*x + parameters(3);
% plot([0:0.01:0.7],f([0:0.01:0.7]),[0;ah_min;0.7],[beta_0;beta_min;beta_07])



end
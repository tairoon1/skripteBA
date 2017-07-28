%%
close all
clc
clear

%% design parameters
% geometry
ll = 100;  %[mm]
hh = ll/8;  %[mm]
bb = hh/2;  %[mm]


% Paris
K_IC = 50; %MPa*m^-0.5
a_init = 0.00001; %[mm]
C = 10^-12; %[?]
m = 2.85; %[?]

for i=1:200
    % load
    P = i; %[N]
    
    %% 1st step: Determine critical crack length
    sigma_max = 3/2*P*ll/bb/hh^2; %[MPa]
    
    beta = @(a) calcBeta(a/hh,1.1,0.15,1.01,1.84); %[-]
    
    a_crit = (K_IC/beta(0)/sigma_max)^2/pi; %[mm]
    
    
    %% 2nd step: Determine fatigue lifex
    
    Nf = paris(2*sigma_max,a_init,a_crit,C,m,beta) %[-]
    % for plotting S-N-curve
    aaa(i)=Nf;
    bbb(i)=sigma_max;
end
% plot
semilogx(aaa,bbb)
xlabel('N')
ylabel('S')
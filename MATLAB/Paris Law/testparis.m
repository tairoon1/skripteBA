Nf=paris(60,0.01,0.03,3.0e-10,3,1.12)

dadN = @(DeltaK) 3e-10*DeltaK.^3

K=0:100;

loglog(K,dadN(K))
grid on
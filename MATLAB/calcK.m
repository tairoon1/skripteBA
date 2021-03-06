Y=@(x) 1.12-0.23.*x+10.6.*x.^2-21.7.*x.^3+30.4.*x.^4;

a = [0.01 0.02 0.03 0.04 0.05];
W = 0.5;
sigma = [96 72 64 60 52];
Y(a/W)
K=Y(a/W).*sigma.*sqrt(pi.*a)
f = @(x) sqrt(x).*Y(x/W);
plot(a,sigma,[0.01:0.001:0.1],mean(K)/sqrt(pi)./f([0.01:0.001:0.1]))
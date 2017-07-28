close all
clc
clear

S56 = 10*[1 330 438 508 558 592 616 636 651];
S58 = 10*[1 298 396 459 504 532 554 571];
S60 = 10*[1 269 357 414 451 475 494];
S62 = 10*[1 244 324 376 408 431 454];
S64 = 10*[1 222 295 341 369 387];
S66 = 10*[1 202 268 308 333 354];


a=0.01:0.005:0.5;

figure(1)
semilogx(S56,a(1:length(S56)))
hold on
semilogx(S58,a(1:length(S58)))
semilogx(S60,a(1:length(S60)))
semilogx(S62,a(1:length(S62)))
semilogx(S64,a(1:length(S64)))
semilogx(S66,a(1:length(S66)))

xlabel('Stress cycles $N$','Interpreter','latex')
ylabel('Crack length $a$','Interpreter','latex')
title('$\log N$-$a$ Plot','Interpreter','latex')
h=legend('$\Delta \sigma=56MPa$','$\Delta \sigma=58MPa$','$\Delta \sigma=60MPa$','$\Delta \sigma=62MPa$','$\Delta \sigma=64MPa$','$\Delta \sigma=66MPa$');
set(h,'Interpreter','latex')
grid on

sample = @(x) 10e-10*x.^3;
beta = @(x) 10.2857*x.^2 - 0.0571*x + 1.1;

for i=1:length(S56)-1
    dadN56(i) = (a(i+1)-a(i))/(S56(i+1)-S56(i));
    dK56(i) = beta(a(i)/0.5)*56*sqrt(pi*a(i));
end

for i=1:length(S58)-1
    dadN58(i) = (a(i+1)-a(i))/(S58(i+1)-S58(i));
    dK58(i) = beta(a(i)/0.5)*58*sqrt(pi*a(i));
end

for i=1:length(S60)-1
    dadN60(i) = (a(i+1)-a(i))/(S60(i+1)-S60(i));
    dK60(i) = beta(a(i)/0.5)*60*sqrt(pi*a(i));
end

for i=1:length(S62)-1
    dadN62(i) = (a(i+1)-a(i))/(S62(i+1)-S62(i));
    dK62(i) = beta(a(i)/0.5)*62*sqrt(pi*a(i));
end

for i=1:length(S64)-1
    dadN64(i) = (a(i+1)-a(i))/(S64(i+1)-S64(i));
    dK64(i) = beta(a(i)/0.5)*64*sqrt(pi*a(i));
end

for i=1:length(S66)-1
    dadN66(i) = (a(i+1)-a(i))/(S66(i+1)-S66(i));
    dK66(i) = beta(a(i)/0.5)*66*sqrt(pi*a(i));
end



% regression
A = [];
b = zeros(length(dadN60)-1,1);
for i=2:length(dadN60)
    b(i-1) = log(dadN60(i));
    A = vertcat(A,[1 log(dK60(i))]);
end
fittedPar=A\b;
fitFun=@(dK) exp(fittedPar(1))*dK.^fittedPar(2);
C=exp(fittedPar(1))
M=fittedPar(2)
% plot DeltaK-dadN
figure(2)
loglog(10:100,sample(10:100));
hold on
loglog(dK56,dadN56)
loglog(dK58,dadN58)
loglog(dK60,dadN60)
loglog(10:100,fitFun(10:100))
loglog(dK62,dadN62)
loglog(dK64,dadN64)
loglog(dK66,dadN66)
xlabel('$\Delta K$','Interpreter','latex')
ylabel('$\frac{da}{dN}$','Interpreter','latex')
title('$\log \Delta K$-$\log\frac{da}{dN}$ Plot','Interpreter','latex')
h=legend('Sample','$\Delta \sigma=56MPa$','$\Delta \sigma=58MPa$','$\Delta \sigma=60MPa$','$Fitted$','$\Delta \sigma=62MPa$','$\Delta \sigma=64MPa$','$\Delta \sigma=66MPa$');
set(h,'Interpreter','latex')
grid on




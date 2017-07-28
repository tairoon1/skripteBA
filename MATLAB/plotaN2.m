close all
clc
clear
N13 =2000*[0 24 39 54 68 81 94 107 119 130 142 154 165];

a=0.12:0.005:0.33;



figure(1)
semilogx(N13,a(1:length(N13)))
hold on
xlabel('Stress cycles $N$','Interpreter','latex')
ylabel('Crack length $a$','Interpreter','latex')
title('$\log N$-$a$ Plot','Interpreter','latex')
h=legend('$\Delta sigma=52MPa$');
set(h,'Interpreter','latex')
grid on

for i=1:length(N13)-1
    dadn13(i) = (a(i+1)-a(i))/(N13(i+1)-N13(i));
    dK13(i) = 1.12*52*sqrt(pi*a(i));
end

figure(2)
loglog(dK13,dadn13)
hold on
xlabel('$\Delta K$','Interpreter','latex')
ylabel('$\frac{da}{dN}$','Interpreter','latex')
title('$\log \Delta K$-$\log\frac{da}{dN}$ Plot','Interpreter','latex')
h=legend('$\Delta sigma=52MPa$');
set(h,'Interpreter','latex')
grid on

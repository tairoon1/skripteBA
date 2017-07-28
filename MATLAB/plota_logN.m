close all
clc
clear
N8 =2000*[0 101 164 224 280 334 387 439 489 539 588 636 683];
N10=2000*[0 52 84 115 144 172 199 226 252 278 303 328 352 376 399 420 441 462 483 504 525 546 566 587 607 626 646 665 684 702 ...
    720 737 754 771 787 803 819 835 850 864 877 890 903];
N14=2000*[0 19 31 42 53 63 73 83 92 101 110 118 126 134 141 149 156 164 172 180 188 195 203 210 217 224 232 239 245 252 259 ...
    265 272 279 285 291 297 303 309 315 321 326 332 338 343 349 354 360 365 370 375 381 386 391 396 401 406 411 416 421 426 ...
    431 436 441 446 451 456 461];

a=0.012:0.0005:0.05;



figure(1)
semilogx(N8,a(1:length(N8)))
hold on
semilogx(N10,a(1:length(N10)))
semilogx(N14,a(1:length(N14)))
xlabel('Stress cycles $N$','Interpreter','latex')
ylabel('Crack length $a$','Interpreter','latex')
title('$\log N$-$a$ Plot','Interpreter','latex')
h=legend('$\Delta sigma=32MPa$','$\Delta sigma=40MPa$','$\Delta sigma=56MPa$');
set(h,'Interpreter','latex')
grid on

for i=1:length(N8)-1
    dadn8(i) = (a(i+1)-a(i))/(N8(i+1)-N8(i));
    dK8(i) = 1.12*32*sqrt(pi*a(i));
end
for i=1:length(N10)-1
    dadn10(i) = (a(i+1)-a(i))/(N10(i+1)-N10(i));
    dK10(i) = 1.12*40*sqrt(pi*a(i));
end
for i=1:length(N14)-1
    dadn14(i) = (a(i+1)-a(i))/(N14(i+1)-N14(i));
    dK14(i) = 1.12*56*sqrt(pi*a(i));
end

figure(2)
loglog(dK8,dadn8)
hold on
loglog(dK10,dadn10)
loglog(dK14,dadn14)
xlabel('$\Delta K$','Interpreter','latex')
ylabel('$\frac{da}{dN}$','Interpreter','latex')
title('$\log \Delta K$-$\log\frac{da}{dN}$ Plot','Interpreter','latex')
h=legend('$\Delta sigma=32MPa$','$\Delta sigma=40MPa$','$\Delta sigma=56MPa$');
set(h,'Interpreter','latex')
grid on

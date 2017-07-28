close all
clc
clear

N11 = 10*[0 79 130 179 227 272 314 356 397 437 477 515 553 592 629 666 702 739 775 810 845 880 914 948 982 1015 1047 1080 1112 1143 1174 1205 1234 1264 1294 1322 1351 1378 1406 1433 1459 1485 1512 1537 1561 1586 1610 1633 1657 1679 1702 1725 1748 1770 1791 1813 1834 1855 1876 1897 1918 1939 1959 1980];
N12 = 10*[0 60 100 138 176 209 242 275 308 337 368 399 429 459 488 517 546 574 602 629 656 682 708 734 759 785 810 835 859 883 906 930 952 975 997 1020 1041 1063 1083 1104 1124 1144 1164 1183 1203 1222 1240 1258 1277 1294 1312 1329 1346 1363 1379 1396 1412 1429 1445 1461 1476 1492 1509 1525 1540];
N13_1 = [0 479 789 1096 1380 1657 1908 2151 2388 2632 2871 3105 3340 3575 3812 4030 4249 4468 4676 4893 5102 5308 5516 5719 5918 6116 6308 6503 6694 6879 7066 7245 7427 7604 7773 7947 8108 8275 8442 8599 8761 8911 9062 9217 9364 9509 9652 9795 9932 10072 10203 10337 10471 10601 10734 10863];
N13 = 10*[0 47 78 109 138 165 190 216 242 266 291 314 337 360 382 405 427 449 471 493 514 535 556 576 596 616 636 656 674 693 712 731 749 766 784 802 818 835 851 867 883 899 914 929 944 959];
N14 = 10*[0 38 63 87 110 131 152 173 193 213 232 251 270 289 307 325 342 360 377 394 411 429 445 461 478 494 509 525 540 554 569 584 598 612 626 640 653 667 680 692 705 717 729 741 753 764 776 788 799 809 821 831 842 852 863 873 883 893 903 914 923 933 944 955 964 975 985 995 1006 1018 1028 1039 1051 1063 1075];
N15 = 10*[0 31 51 71 89 107 124 141 157 173 189 204 219 233 248 263 277 291 305 319 333 346 360 372 385 398 411 423 435 447 459 470 482 493 505 516 526 537 547 558 568 578 588 598 607 617 626 636 645 654 662 671 680 689 697 706 714 722 731 739 747 755 764 772 780 788 797 806];


a=0.12:0.005:0.5;

figure(1)
semilogx(N11,a(1:length(N11)))
hold on
semilogx(N12,a(1:length(N12)))
semilogx(N13_1,a(1:length(N13_1)))
semilogx(N13,a(1:length(N13)))
semilogx(N14,a(1:length(N14)))
semilogx(N15,a(1:length(N15)))

xlabel('Stress cycles $N$','Interpreter','latex')
ylabel('Crack length $a$','Interpreter','latex')
title('$\log N$-$a$ Plot','Interpreter','latex')
h=legend('$\Delta \sigma=44MPa$','$\Delta \sigma=48MPa$','$\Delta \sigma=52MPa$ f=1','$\Delta \sigma=52MPa$','$\Delta \sigma=56MPa$','$\Delta \sigma=60MPa$');
set(h,'Interpreter','latex')
grid on

muster = @(x) 10e-10*x.^3;
beta = @(x) 10.2857*x.^2 - 0.0571*x + 1.1;

for i=1:length(N11)-1
    dadn11(i) = (a(i+1)-a(i))/(N11(i+1)-N11(i));
    dK11(i) = beta(a(i)/0.5)*44*sqrt(pi*a(i));
end

for i=1:length(N12)-1
    dadn12(i) = (a(i+1)-a(i))/(N12(i+1)-N12(i));
    dK12(i) = beta(a(i)/0.5)*48*sqrt(pi*a(i));
end

for i=1:length(N13_1)-1
    dadn13_1(i) = (a(i+1)-a(i))/(N13_1(i+1)-N13_1(i));
    dK13_1(i) = beta(a(i)/0.5)*52*sqrt(pi*a(i));
end

for i=1:length(N13)-1
    dadn13(i) = (a(i+1)-a(i))/(N13(i+1)-N13(i));
    dK13(i) = beta(a(i)/0.5)*52*sqrt(pi*a(i));
end

for i=1:length(N14)-1
    dadn14(i) = (a(i+1)-a(i))/(N14(i+1)-N14(i));
    dK14(i) = beta(a(i)/0.5)*56*sqrt(pi*a(i));
end

for i=1:length(N15)-1
    dadn15(i) = (a(i+1)-a(i))/(N15(i+1)-N15(i));
    dK15(i) = beta(a(i)/0.5)*60*sqrt(pi*a(i));
end


% regression
A = [];
b = zeros(length(dK13),1);
for i=1:length(dK13)
    b(i) = log(dadn13(i));
    A = vertcat(A,[1 log(dK13(i))]);
end
fittedPar=A\b;
fitFun=@(dK) exp(fittedPar(1))*dK.^fittedPar(2);

% plot DeltaK-dadN
figure(2)
loglog(10:1000,muster(10:1000));
hold on
loglog(dK11,dadn11)
loglog(dK12,dadn12)
loglog(dK13_1,dadn13_1)
loglog(dK13,dadn13)
loglog(10:1000,fitFun(10:1000))
loglog(dK14,dadn14)
loglog(dK15,dadn15)
xlabel('$\Delta K$','Interpreter','latex')
ylabel('$\frac{da}{dN}$','Interpreter','latex')
title('$\log \Delta K$-$\log\frac{da}{dN}$ Plot','Interpreter','latex')
h=legend('Paris Law with $C=10\times10^{-10},M=3$','$\Delta \sigma=44MPa$','$\Delta \sigma=48MPa$','$\Delta \sigma=52MPa$ f=1','$\Delta \sigma=52MPa$','$\Delta \sigma=52MPa$ fit','$\Delta \sigma=56MPa$','$\Delta \sigma=60MPa$');
set(h,'Interpreter','latex')
grid on




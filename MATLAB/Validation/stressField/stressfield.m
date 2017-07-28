close all
clc
clear

sigma_inf = 56; %[MPa]

sigmay = @(x,y) 1.12*sigma_inf*sqrt(pi*0.012)/sqrt(2*pi*x)*sqrt(cos(atan(y/x)))*...
                cos(atan(y/x)/2)*(1+sin(atan(y/x)/2)*sin(3*atan(y/x)/2));

[X,Y]=meshgrid(0:0.0005:0.05,-0.025:0.0005:0.025);
values = zeros(size(X,1),size(X,2));

figure
hold on
x=[];
y=[];
for i=1:size(X,1)
    for j=1:size(X,2)
        if X(i,j)<=0.012
            continue
        end
        values(i,j) = sigmay(X(i,j)-0.012,Y(i,j));
        if values(i,j)>=90
            x = [x X(i,j)];
            y = [y Y(i,j)];
        end
    end
end
scatter(x-0.012,y,'rx')
axis([0,0.05,-0.025,0.025])
[x,y,~]=textread('positionStressField.txt');
scatter(x+0.0135,y,'bx')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('Comparison Stress Field','Interpreter','latex')
h=legend('LEFM','LAMMPS');
set(h,'Interpreter','latex')
saveas(gcf,'stressField','epsc')
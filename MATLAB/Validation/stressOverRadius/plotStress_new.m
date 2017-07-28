[r,stress,~]=textread('stress_new.txt');

r = r+0.135;
stress = stress/1.25e-7/1e6;
f = @(r) 15.36320./sqrt(r);

figure
hold on
plot(r,f(r),'rx')
plot(r,stress,'bx')
xlabel('$r$','Interpreter','latex')
ylabel('$\sigma$','Interpreter','latex')
title('Comparison Stress over Radius','Interpreter','latex')
h=legend('LEFM','LAMMPS');
set(h,'Interpreter','latex')
saveas(gcf,'stressOverRadius_new','epsc')
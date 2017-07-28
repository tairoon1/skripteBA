[r,stress,~]=textread('stress.txt');

r = r+0.0135;
stress = stress/1.25e-10/1e6;
f = @(r) 4.85827./sqrt(r);

figure
hold on
plot(r,f(r),'rx')
plot(r,stress,'bx')
xlabel('$r$','Interpreter','latex')
ylabel('$\sigma$','Interpreter','latex')
title('Comparison Stress over Radius','Interpreter','latex')
h=legend('LEFM','LAMMPS');
set(h,'Interpreter','latex')
saveas(gcf,'stressOverRadius','epsc')
function x = plots(par,x,u_bar,Schrodinger)

if (par.Na == 1 && par.Nb==1) || (par.alpha_under==0 && par.beta_under==0)
[~,x] = simulate_base(par,u_bar,Schrodinger,1,1);
figure()
subplot(2,1,1)
plot(linspace(0,par.T,par.K-1),u_bar','LineWidth',3)
xlabel('Time','Interpreter','latex')
ylabel('Control','Interpreter','latex')
axis('tight')
set(gca,'Fontsize',30)
subplot(2,1,2)
plot(linspace(0,par.T,par.K),x','LineWidth',3)
xlabel('Time','Interpreter','latex')
ylabel('State','Interpreter','latex')
axis('tight')
set(gca,'Fontsize',30)
set(gcf,'units','points','position',[1,84,775,782])
end

if par.Na ~= 1 && par.Nb ~=1 && par.alpha_under~=0 && par.beta_under~=0 
a = linspace(par.alpha_min,par.alpha_max,3);
b = linspace(par.beta_min,par.beta_max,3);
[~,x] = simulate_base(par,u_bar,Schrodinger,a,b);
figure()
subplot(3,1,1)
plot(linspace(0,par.T,par.K-1),u_bar(:,:)','LineWidth',3)
xlabel('Time','Interpreter','latex')
ylabel('Control','Interpreter','latex')
axis('tight')
set(gca,'Fontsize',30)
subplot(3,1,2)
plot(linspace(0,par.T,par.K),x','LineWidth',3) 
ylabel('State','Interpreter','latex')
xlabel('Time','Interpreter','latex')
axis('tight')
set(gca,'Fontsize',30)
a = linspace(par.alpha_min,par.alpha_max,21);
b = linspace(par.beta_min,par.beta_max,21);
[error,x] = simulate_base(par,u_bar,Schrodinger,a,b);
[alpha,beta] = meshgrid(a,b);
ax1=subplot(3,1,3);
[c, h]=contourf(alpha,beta,flipud(error'),'LineWidth',3,'ShowText','on');
clabel(c,h,'FontSize',22)
colormap(ax1,sky)
cb = colorbar();
title(cb,'Error $\qquad \qquad$','Interpreter','latex')%,'Rotation',270)
xlabel('$\alpha$','Interpreter','latex')
ylabel('$\beta$','Interpreter','latex')
xticks([par.alpha_min (par.alpha_min+par.alpha_max)/2 par.alpha_max])
yticks([par.beta_min (par.beta_min+par.beta_max)/2 par.beta_max])
set(gcf,'units','points','position',[217,83,919,782])
set(gca,'Fontsize',30)
end

if  (par.Na==1 && par.Nb~=1) || (par.alpha_under == 0 && par.beta_under ~= 0)
a = 1;
b = linspace(par.beta_min,par.beta_max,3);
[~,x] = simulate_base(par,u_bar,Schrodinger,a,b);
figure()
subplot(3,1,1)
plot(linspace(0,par.T,par.K-1),u_bar','LineWidth',3)
xlabel('Time','Interpreter','latex')
ylabel('Control','Interpreter','latex')
axis('tight')
set(gca,'Fontsize',30)
subplot(3,1,2)
plot(linspace(0,par.T,par.K),x','LineWidth',3)
axis('tight')
xlabel('Time','Interpreter','latex')
ylabel('State','Interpreter','latex')
set(gca,'Fontsize',30)
subplot(3,1,3)
b = linspace(par.beta_min,par.beta_max,21);
[error,x] = simulate_base(par,u_bar,Schrodinger,a,b);
plot(b,error,'LineWidth',3)
xlabel('$\beta$','Interpreter','latex')
ylabel('Error','Interpreter','latex')
axis('tight')
set(gcf,'units','points','position',[217,83,919,782])
set(gca,'Fontsize',30)
end

if (par.Na~=1 && par.Nb==1) || (par.alpha_under ~= 0 && par.beta_under == 0)
a = linspace(par.alpha_min,par.alpha_max,3);
b = 1;
[~,x] = simulate_base(par,u_bar,Schrodinger,a,b);
figure()
subplot(3,1,1)
plot(linspace(0,par.T,par.K-1),u_bar','LineWidth',3)
xlabel('Time','Interpreter','latex')
ylabel('Control','Interpreter','latex')
axis('tight')
set(gca,'Fontsize',30)
subplot(3,1,2)
plot(linspace(0,par.T,par.K),x','LineWidth',3)
ylabel('State','Interpreter','latex')
xlabel('Time','Interpreter','latex')
axis('tight')
set(gca,'Fontsize',30)
subplot(3,1,3)
a = linspace(par.alpha_min,par.alpha_max,21);
[error,x] = simulate_base(par,u_bar,Schrodinger,a,b);
plot(a,error,'LineWidth',3)
ylabel('Error','Interpreter','latex')
xlabel('$\alpha$','Interpreter','latex')
axis('tight')
set(gcf,'units','points','position',[217,83,919,782])
set(gca,'Fontsize',30)
end


end


x=linspace(500/400/2,500-500/400/2,400);
subplot(3,1,1)
load('s2_geoxim1.mat')
plot(x,0.3-State.poro)
hold on
load('s2_geoxim4.mat')
plot(x,0.3-State.poro)
load('s2_geoxim5.mat')
plot(x,0.3-State.poro)
xt=importdata('ERICPHI.txt');
plot(xt(:,1),0.3-xt(:,2)./100,'mo','markersize',3)
xlabel('\itx \rm[m]')
ylabel('\it \phi_{HYD}\rm [-]')
legend('t=20 day','t=1 year', 't=5 year','t=5 year (Geoxim)','edgecolor','none')
%%
subplot(3,1,2)
load('s2_geoxim1.mat')
plot(x,State.s(:,1))
hold on
load('s2_geoxim4.mat')
plot(x,State.s(:,1))
load('s2_geoxim5.mat')
plot(x,State.s(:,1))
xt=importdata('ERICsat.txt');
plot(xt(:,1),xt(:,2)./100,'mo','markersize',3)
ylim([0.0 0.25])
xlabel('\itx \rm[m]')
ylabel('\it S_l\rm [-]')
%%
subplot(3,1,3)
load('s2_geoxim1.mat')
plot(x,State.Tk)
hold on
load('s2_geoxim4.mat')
plot(x,State.Tk)
load('s2_geoxim5.mat')
plot(x,State.Tk)
xt=importdata('ERICT.txt');
plot(xt(:,1),xt(:,2)+273.15,'mo','markersize',3)
xlabel('\itx \rm[m]')
ylabel('\it T\rm [K]')
ylim([-23.15+273.15 47+273])
set(gcf,'position',[100 100 500 600])
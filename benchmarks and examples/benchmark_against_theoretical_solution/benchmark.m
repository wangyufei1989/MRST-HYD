x=linspace(500/400/2,500-500/400/2,400);
p=linspace(5.027e6,3.0316e6,400);
T=zeros(400,1);
rhog=100.3993;
cpg=1206.22;
kappa=1e-13;
kr=0.8;
mu=1e-5;
phi=0.3;
rhor=2650;
cpr=800;
tinj=278.9;% the temperature at the first cell
pinj=5.027e6;% the pressure at the first cell
t0=308.17; %note the initial pressure is not the 308.15 as we set in the beginning, it is 308.17 after the first time step;
mujt=1.5e-5;
swi=0.2;
rhow=1002.4;
cpw=4476.29;
vg=kappa*kr/mu*(5e6-3e6)/(500-500/400);
Vg=rhog*cpg*vg/(phi*(swi*rhow*cpw+(1-swi)*rhog*cpg)+(1-phi)*rhor*cpr);
%%
load('ttheo3.mat')
load('stheo3.mat')
T(x<Vg*t)=tinj-mujt.*(pinj-p(x<Vg*t));
T(x>Vg*t)=t0-mujt.*(Vg*t/(500-500/400)*(5e6-3e6));
plot(x,State.Tk)
hold on
plot(x,T,'--')
t=3600*24*365;
load('ttheo4.mat')
load('stheo4.mat')
T(x<Vg*t)=tinj-mujt.*(pinj-p(x<Vg*t));
T(x>Vg*t)=t0-mujt.*(Vg*t/(500-500/400)*(5e6-3e6));
plot(x,State.Tk)
plot(x,T,'--')
load('ttheo5.mat')
t=1.576646579447810e+08;
load('stheo5.mat')
T(x<Vg*t)=tinj-mujt.*(pinj-p(x<Vg*t));
T(x>Vg*t)=t0-mujt.*(Vg*t/(500-500/400)*(5e6-3e6));
plot(x,State.Tk)
plot(x,T,'--')
xlabel('\itx \rm[m]')
ylabel('\it T\rm [K]')
ylim([-43.15+273.15 47+273])
legend('Numerical; t=20 day','Analytical; t=20 day','Numerical; t=1 year','Analytical; t=1 year','Numerical; t=5 year','Analytical; t=5 year','edgecolor','none')
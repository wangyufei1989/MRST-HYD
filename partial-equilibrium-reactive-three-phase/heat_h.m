function [h,u,dh,du]=heat_h(state,drho)

% h: the specific ennthalpy of [water, gas, HYD, rock]
% u: the specific internal energy of [water, gas, HYD, rock]
% dh: the derivitive of h over pressure and temperature
% du: the derivitive of u over pressure and temperature
num=length(state.Tk);

    
 %% field cell; state
h=zeros(num,4);
u=zeros(num,4);
dh=zeros(num,5);
du=zeros(num,5);
href1=0;href2=0;href3=-368622;href4=0;
% href3=0;
tref1=273.15;tref2=273.15; tref3=273.15;tref4=273.15;
pref=3e6;
cp1=state.heat_capacity(1);
cp2=state.heat_capacity(2);
cp3=state.heat_capacity(3);
cp4=state.heat_capacity(4);
%cp1=3950;cp2=1832;cp3=1930;cp4=800;


mutj=state.mutj*cp2;
h(:,1)=href1+cp1.*(state.Tk-tref1);%water
h(:,2)=href2+cp2.*(state.Tk-tref2)-mutj.*(state.pressure(:,2)-pref);%co2
h(:,3)=href3+cp3.*(state.Tk-tref3);%hyd
h(:,4)=href4+cp4.*(state.Tk-tref4);%rock
u(:,1)=h(:,1)-state.pressure(:,1)./state.rho(:,1);
u(:,2)=h(:,2)-state.pressure(:,2)./state.rho(:,2);

 u(:,3)=h(:,3);
 u(:,4)=h(:,4);
%u=h;
%%
dh(:,1)=zeros(num,1);%dhl/dpl
dh(:,2)=-mutj.*ones(num,1);%dhg/dpg
dh(:,3)=cp1.*ones(num,1);%dhl/dt
dh(:,4)=cp2.*ones(num,1);%dhg/dt
dh(:,5)=cp3.*ones(num,1);%dhh/dt
dh(:,6)=cp4.*ones(num,1);% dhr/dt

du(:,1)=dh(:,1)-1./state.rho(:,1)+state.pressure(:,1)./state.rho(:,1).^2.*drho(:,1);
du(:,2)=dh(:,2)-1./state.rho(:,2)+state.pressure(:,2)./state.rho(:,2).^2.*drho(:,2);
du(:,3)=dh(:,3)+state.pressure(:,1)./state.rho(:,1).^2.*drho(:,5);
du(:,4)=dh(:,4)+state.pressure(:,2)./state.rho(:,2).^2.*drho(:,6);
du(:,5)=dh(:,5);
du(:,6)=dh(:,6);

%du=dh;




end

clear
close all





for  i=[0.5 1 2 3 3.35 3.5 3.8 4 4.5]
    figure
    state.pressure(:,1)=linspace(100*i,100*i+100,100).*1e5;
 state.pressure(:,2)=linspace(100*i,100*i+100,100).*1e5;
state.Tk=linspace(276.15,306.15,100)';

state.m_NaCl=0.5.*ones(100,1);
state.species=[4.23879216936390e-08,1.63392933157338e-07,0,4.23879216936390e-08.*0.5.*0.023,4.23879216936390e-08.*0.5.*0.0355,0,4.23879216936390e-08.*0.06];
[ state.rho,state.mu]=rhomu_p_frac_kinetic_h(state);

rho=state.rho;mu=state.mu;

    plot((1:100)*10,rho(:,1),'-')
hold on
 plot((1:100)*10,rho(:,2),'--')
    
end



%%
clear
close all
Rock.poro=0.3.*ones(100,1);
Rock.perm=1e-13.*ones(100,1);
pJ=sqrt(mean(Rock.poro)./mean(Rock.perm)).*sqrt(Rock.perm./Rock.poro);
% residual saturation (s_lr,s_gr); maximum relative permeability
% (k_lm,k_gm);
s_lr=0.2e0;s_gr=0.05;
k_lm=1;k_gm=0.8;
alpha_van=5;%[bar]
m_van=0.8;
%  van Genuchten model or simple model can be used; in the simple model;
%  the relative permeability is simply proportional to the saturation.
van_genuchtten=1;
pcmax=50e5;

Fluid = initFluid_VG('sr' , [ s_lr, s_gr] , ...
 'kwm', [ k_lm, k_gm], 'pc_scale', pJ,'alpham',[alpha_van,m_van],'pcmax',pcmax);



for i=[1 1.5 2 3 3.5 4 4.5]
    figure
    state.pressure(:,1)=linspace(100*i,100*i+100,100).*1e5;
 state.pressure(:,2)=linspace(100*i,100*i+100,100).*1e5;
state.Tk=linspace(276.15,306.15,100)';

state.m_NaCl=0.5.*ones(100,1);
state.species=[4.23879216936390e-08,1.63392933157338e-07,0,4.23879216936390e-08.*0.5.*0.023,4.23879216936390e-08.*0.5.*0.0355,0,8.69509205359146e-10];
PEQ=Fluid.peq(state);

    plot((1:100)*10,state.pressure(:,1),'-')
hold on
 plot((1:100)*10,PEQ,'--')
    
end

%%
close all
figure

plot([1000, 2000,3000,3500,4000,4500],[0 0 0 170 450 770],'m--','linewidth',3)

hold on

plot([1000, 2000,3000,3500,4000,4500],[35 160 250 280 300 310],'c-','linewidth',3)

legend('NBZ','HFZ','edgecolor','none')

xlabel('Ocean floor depth [m]')
ylabel('Thickness [m]')
%%
close all
figure

plot([500 1000, 2000,3000, 3350,3500,3800,4000,4500],[0 0 0 0 0 30 200 310 580],'m--','linewidth',3)

hold on

plot([500 1000, 1500, 2000,3000,3500,4000,4500],[0 0 70 130 220 250 265 280],'c-','linewidth',3)

legend('NBZ','HFZ','edgecolor','none','fontsize',12)

xlabel('Ocean floor depth [m]')
ylabel('Thickness [m]')
%%
close all
for i=[  2 3.8 4.5]
    figure
    state.pressure(:,1)=linspace(100*i,100*i+100,100).*1e5;
 state.pressure(:,2)=linspace(100*i,100*i+100,100).*1e5;
state.Tk=linspace(276.15,306.15,100)';

state.m_NaCl=0.5.*ones(100,1);
state.species=[4.23879216936390e-08,1.63392933157338e-07,0,4.23879216936390e-08.*0.5.*0.023,4.23879216936390e-08.*0.5.*0.0355,0,4.23879216936390e-08.*0.06];
PEQ=Fluid.peq(state);

    plot(linspace(276.15,306.15,100)',ones(100,1)*i*1e7,'-')
hold on
 plot(state.Tk,PEQ,'--')
    
end

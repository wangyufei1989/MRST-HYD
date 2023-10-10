clear;%clear the workspace
close all % close all figures
%% install the mrst_co2
%add the path of the folder of mrst_co2;for windows system, you can use addpath('\..\..'); or simply find the 'startup.m' file in the main folder and run it. 
addpath('../../') 
run startup.m
%% define the grid system
% grid discretization; the numbers of grids in the (x,y,z) direction are (nx,ny,nz);
nx=400;ny=1;nz=1;
dims= [nx ny nz];
% the domain sizes in the (x,y,z) directions are (distance, thickness, depth);
distance=500; thickness=10; depth=50;
domain=[distance thickness  depth]; 
% generate cartesian grid  system; the dimension can be three;
Grid= computeGeometry(cartGrid(dims,domain));
% define the intrinsic permeability
Rock.perm=ones(nx*ny*nz,1).*1e-13;
% assign uniform porosity to all the grid;
Rock.poro = 0.3*ones( nx*ny*nz,1);
%% define the fluid property; the retension curve, relative permeability and equilibrium pressure;
%  van Genuchten model is used; one can also modify the  model if necessary.
%J-fucntion; Leverette function 
pJ=sqrt(mean(Rock.poro)./mean(Rock.perm)).*sqrt(Rock.perm./Rock.poro);
% residual saturation (s_lr,s_gr); maximum relative permeability (k_lm,k_gm);
s_lr=0.2e0;s_gr=0;
% maximum relative permeability
k_lm=1;k_gm=0.8;
% scaling pressure [bar]
alpha_van=5;
% shape factor
m_van=0.8;
%maximum capillary pressure
pcmax=50e5;
Fluid = initFluid_VG('sr' , [ s_lr, s_gr] , ...
 'kwm', [ k_lm, k_gm], 'pc_scale', pJ,'alpham',[alpha_van,m_van],'pcmax',pcmax);
%% define the fluid density
% if density is set simple (rhomu_flag=1), then we use the model
% rho=rho_ref*(1+c_p(p-p_r)-c_t(T-T_r));mu=constant
rhomu_flag=1;
% give the reference [rhor_lquid,rhor_gas]
rhor=[1000 100];
% give [c_p_liquid, c_p_gas, c_t_liquid, c_t_gas]
compress=[0.5e-9,1e-9,0,0];
%give the reference pressure and temperature
pr=[1e5 1e5];tr=[308.15 308.15];
% give the viscosity
mu=[1e-3,1.04e-5];
rhomu=struct('rhomu_flag',rhomu_flag,'rhor',rhor,'compress',compress,'mu',mu,'pr',pr,'tr',tr);
heat_capacity=[4476.29 1206.22 1930.57 800];
mutj=1.5e-5;
%% define the well location and rate
% in this benchmark, there is no injection well; so we set Well=[]; 
% examples of well definition can be found in
% 'one_dimensional_example_of_deep_ocean_GCS' or 'two_dimensional_example_of_deep_ocean_GCS'
Well=[];
%% initialize the fluid system
% set the gravity; if the gravity is  necessary, then use 'gravity reset on';
gravity reset off
% set the initial molality of aqueous NaCl:  m_NaCl;
m_NaCl=0.5e-4;
% set the initial temperature
Tk=273.15+35-0.85*30/35;
%initial capillary pressure (pc), liquid pressure and gas pressure at the top (p_l0,p_g0); the
%pressures increase downward under gravity force, if gravity is on;
pc=29200;
p_l0=50*barsa;
p_g0=50.*barsa+pc;
% reference density of water, CO2, HYD, halite and rock
density_ref=[1000 640 917 2160 2650];
%dispersion parameters of two phases: [dispersivity_longitudinal,
%dispersivity_transverse,molecular_diffusion];
D_liquid=[0 0 0];D_gas=[0 0 0];
%% thermal conductivity J/(s m K) of liquid gas rock and hydrate
% for fluid phase we may need dispersivity
DT_liquid=[0 0 0];DT_gas=[0 0 0]; 
% for rock and HYD phases
DT_R=[0];DT_HYD=[0];
%initial hydrate fraction and halite fraction in the solid phase; in the current version
%this term is not used.
ss0= [0 0];
% define if we consider the dissolution of CO2; disco2=0 means the
% dissoltion is negligibly small
disco2=0;
pressure_distribution.val=[50*barsa,30*barsa];
pressure_distribution.type='linear';
State=initStatePP_kinetic_unstructure_ocean(Grid,Well,[p_l0,p_g0],density_ref,Tk,Fluid,m_NaCl,Rock,D_liquid,D_gas,...
    DT_liquid,DT_gas,DT_R,DT_HYD,ss0,rhomu,disco2,heat_capacity,mutj,pressure_distribution); 
% if you want to use the initial data from observation, modify the 'State' at
% here; e.g., State.pressure(:,1)=[...], State.pressure(:,2)=[...],

%%  pressure boundary
% the boundnary pressure now only allow the pressure boundary.
% if bc=[], then all the boundaries are closed. here 'xmax' means the right
% most boundary; similarly,'zmax' means the bottom boundary; we give the
% pressure at the right side;

% add the right boundary with liquid pressure and gas pressure (p_l,p_g)
pB_out=[30*barsa,30.*barsa+pc];
% the mass fraction of water component, CO2 component, Cl component and Z
% component in liquid phase on the boundary; 
frac_out=[1 1 1 1];
% the boundary temperature
T_out=308.15-0.75;
% the relative permeability on the boundary if flow goes in side of the
% domain;
kr=[1 1];
bc  = pside_h([], Grid, 'xmax', pB_out,frac_out,T_out,kr,State,Fluid); 

%add the left boundary (the same as we add the right boundary)
pB_in=[30*barsa,30*barsa+pc]+20*barsa*1;
frac_in=[1 1 1 1 ];
T_in=278.15-0.75;
kr=[0 1];
bc= pside_h(bc, Grid, 'xmin', pB_in,frac_in,T_in,kr,State,Fluid); 
bc.species=[1 1 0 1*m_NaCl*0.023 1*m_NaCl*0.0355 0 0];
omega_t=[1,0,0,0,0,0.710526315789474,0;0,1,0,0,0,0.289473684210526,1;0,0,1,0,1.64788732394366,0,0;0,0,0,1,-0.647887323943662,0,0];
bc.omega_t=omega_t;
bc.frac=(omega_t(:,[1 4 5 7])*bc.species([1 4 5 7])'./sum(bc.species([1 4 5 7])))';
bc.frac=repmat(bc.frac,size(bc.value,1),1);
bc.salinity=sum(bc.species([ 4 5 7]))/sum(bc.species([1 4 5 7]));
bc.m_NaCl=0.5.*(bc.species(4)/0.023+bc.species(5)/0.0355)/bc.species(1);
bc.heat_capacity=heat_capacity;
bc.mutj=mutj;
for i=1:length(bc.face)
     bc.type{i}='pressure';
end
for i=1:length(bc.face)
    
    if bc.type{i}=='pressure'
        bc.tb(i)=1;
    else
        bc.tb(i)=0;
    end
    
end
%% kinetic rate law 
%r_k=kinetic_rate*specific_area*reactive_fraction*exp(-Ea/RT)(p_g-p_eq)
% intrinsic kinetic rate
kinetic_rate=5.4e-13;
% activation energy for hydrate formation and dissociation
Ea1=0;Ea2=0;
State.kinetic_rate=kinetic_rate;
State.kinetic_rate0=kinetic_rate;
State.Ea1=Ea1;
State.Ea2=Ea2;
% if the molar abundance is smaller than nr_chemical the reaction rate is
% zero.
State.nr_chemical=1000*1e-3.*Grid.cells.volumes.*0.1./0.018;

State.permr=Rock.perm(:,1);
% we do not use Kozeny Carman model
State.kozeny_carman=0;

%% simulation time and initial time step 
% assign the total simulation time 
%T= distance/darcyv*1e6;
T=3600*24*365*5.1;
% initial time step.
dt=1e-2;
% maximum allowed time step
State.dtmax=1E4;
%record the time step;
DT=[];
% store frequency 
dst=3600*24*10;
% indicator for storage times
jc=1;
% initialize time
t=0;
%% error monitor
% tolerance is based on the maximum allowed residual pressure; when the
% pressure change is smaller than dplim, we deem the result converges.
State.dplim=1e-3;
State.pcmax=pcmax;
% relax the tolerance if necessary, 1.0 means there is not relaxation; 10.0
% means dplim is increased by 10 times
State.dplim_relax=50;
% indicator for convergence solve=0 means does not converge
State.solve=1;
% initial mass of the CO2 component
mc=sum((State.rho(:,1).*State.s(:,1).*State.frac(:,2)).*State.poro.*Grid.cells.volumes)+sum((State.rho(:,2).*State.s(:,2)).*State.poro.*Grid.cells.volumes...
    +sum(State.rhoref(:,3).*(0.3-State.poro).*Grid.cells.volumes).*44/(44+18*6));
% initial mass of water component
mb=sum((State.rho(:,1).*State.s(:,1).*State.frac(:,1)).*State.poro.*Grid.cells.volumes)+sum(State.rhoref(:,3).*(0.3-State.poro).*Grid.cells.volumes).*18*6/(44+18*6);
% total CO2 influx
co2=0;
% total water influx
water=0;
%% main loop

while t<T
    
    
    % the Newton loop for update all state variables
    [State,Well,dt]  =  Newton_loop(State, Grid, Fluid,Rock ,dt,'wells', Well,'bc',bc);

    %update time if the simulation converges
    if State.solve==1
        t=t+dt;
    end
    % save results
    if  abs(t-dst)<3*dt && State.solve==1&&t<T
        dissolve=(State.poro-State.poro0)./dt;
        savei=num2str(jc);
        savex=strcat('s21_geoxim',savei,'.mat');
        save(savex,'State')
        savet=strcat('t21_geoxim',savei,'.mat');
        save(savet,'t')
        jc=jc+1;
        
        if jc==2
            dst=3600*24*20;
        end
        
        if jc==3
            dst=3600*24*180;
        end
        
        if jc==4
            dst=3600*24*365;
        end
        
        if jc==5
            dst=3600*24*365*5;
        end
        
        
    end
    % next out flux of CO2 component through the boundary
    co2=co2-State.outflux(2);
    % next out flux of water component through the boundary
    water=water-State.outflux(1);
    % accumulative error for the CO2 component; mass balance for the CO2 component
    err2=(sum((State.rho(:,1).*State.s(:,1).*State.frac(:,2)).*State.poro.*Grid.cells.volumes)+sum((State.rho(:,2).*State.s(:,2)).*State.poro.*Grid.cells.volumes)...
        +sum(State.rhoref(:,3).*(0.3-State.poro).*Grid.cells.volumes).*44/(44+18*6)-co2 -mc)...
        /(sum((State.rho(:,1).*State.s(:,1).*State.frac(:,2)).*State.poro.*Grid.cells.volumes)+sum((State.rho(:,2).*State.s(:,2)).*State.poro.*Grid.cells.volumes)...
        +sum(State.rhoref(:,3).*(0.3-State.poro).*Grid.cells.volumes).*44/(44+18*6));
    
    % accumulative error for the water component; mass balance for the water component
    err1=(sum((State.rho(:,1).*State.s(:,1).*State.frac(:,1)).*State.poro.*Grid.cells.volumes)+sum(State.rhoref(:,3).*(0.3-State.poro).*Grid.cells.volumes).*18*6/(44+18*6)...
        -water -mb)...
        /(sum((State.rho(:,1).*State.s(:,1).*State.frac(:,1)).*State.poro.*Grid.cells.volumes)+sum(State.rhoref(:,3).*(0.3-State.poro).*Grid.cells.volumes).*18*6/(44+18*6));

    %% update the time step
    DT=[DT;dt];
    State.meandt=mean(DT);
    dt=Updatedt_kinetic(State,dt)
    %% plot result
    q=1;
    if q==1
    subplot(3,1,1)
    plot(Grid.cells.centroids(1:nx),(State.poro(:,1)))
    xlabel('x (m)')
    ylabel('porosity')    
    subplot(3,1,2)
    plot(Grid.cells.centroids(1:nx),State.s(:,2))
    xlabel('x (m)')
    ylabel('gas saturation')
    subplot(3,1,3)
    plot(Grid.cells.centroids(1:nx),State.Tk(:))
    xlabel('x (m)')
    ylabel('temperature')
    drawnow
    end
    %% only for debug
    if isnan(sum(State.peq(:)))||isnan(sum(State.species(:)))
        pause
    end
end


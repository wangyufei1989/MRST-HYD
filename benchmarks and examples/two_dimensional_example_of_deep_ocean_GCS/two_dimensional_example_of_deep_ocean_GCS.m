clear;%clear the workspace
close all % close all figures
%% install the mrst_co2
%add the path of the folder of mrst_co2;for windows system, you can use addpath('\..\..'); or simply find the 'startup.m' file in the main folder and run it. 
addpath('../../') 
run startup.m
%% define the grid system
% grid discretization; the numbers of grids in the (x,y,z) direction are (nx,ny,nz);
nx=20;ny=1;nz=40;
dims= [nx ny nz];
% the domain sizes in the (x,y,z) directions are (distance, thickness, depth);
distance=1000;thickness=60; depth=400;
domain=[distance thickness  depth]; 
% generate cartesian grid  system; the dimension can be three;
Grid= computeGeometry(cartGrid(dims,domain));
% change the caresian grid system to radial system
Grid=orth2radial(Grid);
% horizontal to vertical permeability ratio
ani_ratio=5; 
% vertical permeability [m^2]
kv=1e-12;
% assign the vertical and horizontal permeability [kx,ky;kz]
Rock.perm=[ones(nx*ny*nz,1).*kv.*ani_ratio,ones(nx*ny*nz,1).*kv.*ani_ratio,ones(nx*ny*nz,1).*kv];
% assign uniform porosity to all the grid;
Rock.poro = 0.3*ones( nx*ny*nz,1);
Rock.ani=ani_ratio;

%% define the fluid property; the retension curve, relative permeability and equilibrium pressure;
%  van Genuchten model is used; one can also modify the  model if necessary.
%J-fucntion; Leverette function 
pJ=sqrt(mean(Rock.poro)./mean(Rock.perm(:,3))).*sqrt(Rock.perm(:,3)./Rock.poro);
% residual saturation (s_lr,s_gr); maximum relative permeability (k_lm,k_gm);
s_lr=0.2;s_gr=0.05; k_lm=0.8;k_gm=1;
% scaling pressure [bar]
alpha_van=5;
% shape factor
m_van=0.8;
%maximum capillary pressure
pcmax=50e5;
% generate the object storing the function for relative permeability,
% retension curve and hydration equilibrium pressure; the equilibrium
% pressure is a function of temperature and salinity based on interpolation
Fluid = initFluid_VG('sr' , [ s_lr, s_gr] , ...
 'kwm', [ k_lm, k_gm], 'pc_scale', pJ,'alpham',[alpha_van,m_van],'pcmax',pcmax);

%% define the fluid density
% we set rhomu_flag=0, if density and viscosity use the sophisticated model given by ''Numerical modeling of geological carbon  
% sequestration: enhanced dissolution in randomly heterogeneous media (https://upcommons.upc.edu/handle/2117/376077)'' 
rhomu_flag=0;
rhomu=struct('rhomu_flag',rhomu_flag);
heat_capacity=[4476.29 1206.22 1930.57 800];
mutj=0;
%% define the well location and rate
% the injected volume rate with the reference density 640 kg/m^3; for example 1 Mt/year
rate=1e6/.64/24/365/3600;
% injected phase; [0 1] means inject gas; [1 0] means inject brine;
phase=[0 1];
% well radius 
r_w=0.1;
% add the vertical injection well;  the (x,y) location of the well is
% (wx,wy), which means the x location of the well is the wx-th grid; the y
% location of the well is the wy-th grid; the well pierces from the wz1
% layer to the wz2 layer in the vertical direaction; 1<=wx<=nx; 1<=wy<=ny; 1<=wz1<=wz2<=nz;
wx=1;wy=1;wz1=nz-3;wz2=nz-0;
% injection well temperature
W_Tk=14+273.15;
% initial well pressure; normally, this value should be slightly larger or equal to the initial pressure; 
pW=40.11*barsa;
% molality of salt in  the water in the injection well
W_m_NaCl=1e-2;
W_species=[1 1 0 W_m_NaCl*0.023 W_m_NaCl*0.0355 0 0];
% the kernal matrix that transfers species to component u_C=omega_t*xi
omega_t=[1,0,0,0,0,0.710526315789474,0;0,1,0,0,0,0.289473684210526,1;0,0,1,0,1.64788732394366,0,0;0,0,0,1,-0.647887323943662,0,0];
% mass fraction of four components (w, c, Cl, Z) in water in the injection well
W_frac=omega_t(:,[1 4 5 7])*W_species([1 4 5 7])'./sum(W_species([1 4 5 7]));
% salinity of the water in the injection well
W_sal=sum(W_species([ 4 5 7]))/sum(W_species([1 4 5 7]));


Well = verticalWell_h([],Grid, Rock,wx,wy, (wz1:wz2), 'Type', 'rate', ...
    'Val', rate,'name', 'i1', 'radius', r_w, 'Comp_i', phase,'pressure',[pW pW],'species',W_species,...
    'frac',W_frac,'Tk',W_Tk,'salinity',W_sal,'m_NaCl',W_m_NaCl);
Well.heat_capacity=heat_capacity;
Well.mutj=mutj;

% if you need more than one injection well, then use
%Well = verticalWell_h(Well,Grid, Rock,wx,wy, (wz1:wz2), 'Type', 'rate', ...
 %   'Val', rate,'name', 'i1', 'radius', r_w, 'Comp_i', phase,'pressure',[pW pW],'species',W_species,...
  %  'frac',W_frac,'Tk',W_Tk,'salinity',W_sal,'m_NaCl',W_m_NaCl);
%% initialize the fluid system
% set the gravity; if the gravity is not necessary, then use 'gravity reset off';
gravity reset on
% set the molality of NaCl:  m_NaCl;
m_NaCl=0.5;
% set the temperature
Tk=273.15+3;
%initial liquid pressure and gas pressure at the top (p_l0,p_g0); the
%pressures increase downward under gravity force;
p_l0=350*barsa;
p_g0=1.*barsa;
% reference density of water, CO2, HYD, halite and rock
density_ref=[1000 640 917 2160 2650];
%dispersion parameters of two phases: [dispersivity_longitudinal,
%dispersivity_transverse,molecular_diffusion];
D_liquid=[5 0.5 1e-9];D_gas=[0 0 1e-9];

%% thermal conductivity J/(s m K) of liquid gas rock and hydrate
% thermal conductivity of water and CO2; for fluid phase we may need dispersivity
DT_liquid=[0 0 0.6];DT_gas=[0 0 0.14]; 
% thermal conductivity of rock and hydrate
DT_R=[1.5];DT_HYD=[0.5];
%initial hydrate fraction and halite fraction in the solid phase; in the current version
%this term is not used.
ss0= [0 0];
% define if we consider the dissolution of CO2; disco2=0 means the
% dissoltion is negligibly small
disco2=1;
State=initStatePP_kinetic_unstructure_ocean(Grid,Well,[p_l0,p_g0],density_ref,Tk,Fluid,m_NaCl,Rock,D_liquid,D_gas,...
    DT_liquid,DT_gas,DT_R,DT_HYD,ss0,rhomu,disco2,heat_capacity,mutj); 
% if you want to use the initial data from observation, modify the 'State' at
% here; e.g., State.pressure(:,1)=[...], State.pressure(:,2)=[...],

%%  pressure boundary
% the boundnary pressure now only allow the pressure boundary.
% if bc=[], then all the boundaries are closed. here 'xmax' means the right
% most boundary; similarly,'zmax' means the bottom boundary; we give the
% pressure at the right side;

% gas pressure and liquid pressure at the top 
pB_top=[350*barsa,0.99*barsa];
% we consider that the four components in water
frac=[1 1 1 1];
T_top=Tk;
% the boundary is permeable to both water and CO2
kr=[1 1];

% out boundary
bc  = pside_h([], Grid, 'xmax', pB_top,frac,T_top,kr,State,Fluid); 

% top boundary
bc= pside_h(bc, Grid, 'zmin', pB_top,frac,T_top,kr,State,Fluid); 
  
 
  
% the gas pressure on the boundary is constant
bc.value(:,2)=ones(size(bc.value(:,2))).*pB_top(2);
% we only need the mass fractions of four components and salinity in the liquid phase
bc.species=[1 1 0 1*m_NaCl*0.023 1*m_NaCl*0.0355 0 0];
bc.omega_t=[1,0,0,0,0,0.710526315789474,0;0,1,0,0,0,0.289473684210526,1;0,0,1,0,1.64788732394366,0,0;0,0,0,1,-0.647887323943662,0,0];
bc.frac=(omega_t(:,[1 4 5 7])*bc.species([1 4 5 7])'./sum(bc.species([1 4 5 7])))';
bc.frac=repmat(bc.frac,size(bc.value,1),1);
bc.salinity=sum(bc.species([ 4 5 7]))/sum(bc.species([1 4 5 7]));
bc.m_NaCl=0.5.*(bc.species(4)/0.023+bc.species(5)/0.0355)/bc.species(1);
bc.heat_capacity=heat_capacity;
bc.mutj=mutj;
% indicator for constant pressure boundary
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
kinetic_rate=1e-13;
% activation energy for hydrate formation and dissociation
Ea1=0;Ea2=0;
State.kinetic_rate=kinetic_rate;
State.kinetic_rate0=kinetic_rate;
State.Ea1=Ea1;
State.Ea2=Ea2;
% if the molar abundance is smaller than nr_chemical the reaction rate is
% zero.
State.nr_chemical=1000*1e-3.*Grid.cells.volumes.*0.1./0.018;
State.permr=Rock.perm(:,3);
% we do not use Kozeny Carman model
State.kozeny_carman=0;

%% simulation time and initial time step 
% assign the total simulation time 
T=3600*24*365*20;
% assign the injection time
T_INJ=3600*24*365*20;
% initial time step. very small. 
dt=1e-2;
% maximum allowed time step
State.dtmax=1E4;
%record the time step;
DT=[];
% indicator for storage times
jc=1;
% initialize time
t=0;
% store frequency 
dst=3600*24*365*1;
%% error monitor
% tolerance is based on the maximum allowed residual pressure; when the
% pressure change is smaller than dplim, we deem the result converges.
State.dplim=1e-3;
% relax the tolerance if necessary, 1.0 means there is not relaxation; 10.0
% means dplim is increased by 10 times
State.dplim_relax=50;
State.pcmax=pcmax;
% indicator for convergence solve=0 means does not converge
State.solve=1;
% initial mass of CO2 component
mc=sum((State.rho(:,1).*State.s(:,1).*State.frac(:,2)).*State.poro.*Grid.cells.volumes)+sum((State.rho(:,2).*State.s(:,2)).*State.poro.*Grid.cells.volumes...
    +sum(State.rhoref(:,3).*(0.3-State.poro).*Grid.cells.volumes).*44/(44+18*6));
% initial mass of water component
mb=sum((State.rho(:,1).*State.s(:,1).*State.frac(:,1)).*State.poro.*Grid.cells.volumes)+sum(State.rhoref(:,3).*(0.3-State.poro).*Grid.cells.volumes).*18*6/(44+18*6);
% net mass influx of CO2 component 
co2=0;
% net mass influx of water component
water=0;
% injected CO2 mass
inj=0;
err=[];
i=1;
while t<T 
    i=i+1;
    %% loop for the reactive transport; update all state variables
    [State,Well,dt]  =   Newton_loop(State, Grid, Fluid,Rock ,dt,'wells', Well,'bc',bc);
    if State.solve==1
        t=t+dt;     
    end
    %% monitor the accumulated mass balance of water component and CO2 component
    inj=inj+Well.val*640*dt;
    co2=co2+Well.val*640*dt-State.outflux(2);
    water=water-State.outflux(1);
    err2=(sum((State.rho(:,1).*State.s(:,1).*State.frac(:,2)).*State.poro.*Grid.cells.volumes)+sum((State.rho(:,2).*State.s(:,2)).*State.poro.*Grid.cells.volumes)...
        +sum(State.rhoref(:,3).*(0.3-State.poro).*Grid.cells.volumes).*44/(44+18*6)-co2 -mc)...
        /(sum((State.rho(:,1).*State.s(:,1).*State.frac(:,2)).*State.poro.*Grid.cells.volumes)+sum((State.rho(:,2).*State.s(:,2)).*State.poro.*Grid.cells.volumes)...
        +sum(State.rhoref(:,3).*(0.3-State.poro).*Grid.cells.volumes).*44/(44+18*6))
    
    
    err1=(sum((State.rho(:,1).*State.s(:,1).*State.frac(:,1)).*State.poro.*Grid.cells.volumes)+sum(State.rhoref(:,3).*(0.3-State.poro).*Grid.cells.volumes).*18*6/(44+18*6)...
        -water -mb)...
        /(sum((State.rho(:,1).*State.s(:,1).*State.frac(:,1)).*State.poro.*Grid.cells.volumes)+sum(State.rhoref(:,3).*(0.3-State.poro).*Grid.cells.volumes).*18*6/(44+18*6))
    State.inj=inj;
    %% save results    
    if  abs(t-dst)<6*dt && State.solve==1
        savei=num2str(jc);
        savex=strcat('s_2dsss','t',savei,'.mat');
        parsave(savex,State,t,Well)
        jc=jc+1;
        dst=jc*year;
    end
   %% update time step length      
    DT=[DT;dt];
    State.meandt=mean(DT);
    [dt,State]=Updatedt_kinetic_ocean(State,dt);
    dt
    
    %% figure
    if fix(i/10)==i/10
        X=Grid.cells.centroids(1:nx);
        Y=linspace(0,depth,nz);
        [X,Y]=meshgrid(X,Y);
        S=reshape(State.s(:,2),nx,nz);
        contourf(X,Y,S', 'LineStyle', 'none')
        xlim([0 distance])
        xlabel('r (m)')
        ylabel('z (m)')
        colorbar
        set (gca,'Ydir','reverse')
        drawnow
    end
    %% for debug
    if isnan(sum(State.peq(:)))||isnan(sum(State.species(:)))
        pause
    end
end
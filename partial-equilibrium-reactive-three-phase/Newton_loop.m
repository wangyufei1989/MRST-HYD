function [state,W,dt] =Newton_loop(state, G, fluid, rock,dt,varargin)
%  main Newton-Raphson loop for PARTIAL equilibrium reactive multiphase phase flow
%  model
%
% SYNOPSIS:
%[state,W,dt] =Newton_loop(state, G, fluid, rock,dt,varargin)
%
% DESCRIPTION:
%   This function serves to update the independent transport variables: lquid
%   pressure (pl), gas pressure (pg), porosity (phi) and temperature (T); and
%   independent reaction variables: water component, CO2 component, Cl component 
%   and charge (Z) component; and the dependent variables.
% REQUIRED PARAMETERS:
%    state    - an objective containing all state variables
%    G        - grid system
%    fluid    - containing relative permeability curves, retention curve and
%    equilibrium hydration pressure curve
%    rock     - containing initial porosity and intrinsic permeability
%    dt       - time step
%    well     - injection well
%    bc       - boundary
% RETURNS:
%    state    - an objective containing all state variables
%    W        - an objective containing information of the injection well
%    dt       - time step may be reduced to a very small value if iteration fails to converge.
% SEE ALSO:
%  kinetic_reaction_hc

%{

This file is part of mrst_ch based on MRST.

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

opt = struct('bc', [], 'src', [], 'wells', [], ...
    'bcp',[],...
    'LinSolve',     @mldivide,        ...
    'MatrixOutput', false, ...
    'Verbose',      mrstVerbose,...
    'gravity',      gravity(), ...
    'condition_number',false,...
    'pc_form','nonwetting',...
    'use_trans',false);

opt = merge_options(opt, varargin{:});
g_vec   = opt.gravity;
grav = (norm(g_vec(1:G.griddim)) > 0) || isfield(G, 'grav_pressure');
if all([~opt.MatrixOutput , ...
        isempty(opt.bc)   , ...
        isempty(opt.src)  , ...
        isempty(opt.bcp)  , ...
        isempty(opt.wells), ~grav])
    warning(msgid('DrivingForce:Missing'),                   ...
        ['No external driving forces present in model--', ...
        'state remains unchanged.\n']);
end
%% use the variable updated from last time step as the initial variable of this new timestep
ncc=G.cells.num;% number of cells.
[state.perm,~]=kdk(rock,state);% updatate the intrinsic permeability based on Kozeny-Carmen model (or other models)
hT=computeTrans(G, state);% update the half transmissivity of each face with respect to each cell
fht=1./accumarray(G.cells.faces(:,1), 1./hT);% calculate the transmissivity of each internel face
% The detailed description of calculation of the 'hT' and the 'fht' can be found in pages 131 -134 of LIE, Knut-Andreas. An introduction to
% reservoir simulation using MATLAB/GNU Octave: User guide for the MATLAB Reservoir Simulation Toolbox (MRST).
% Cambridge University Press, 2019.
nw=length(opt.wells);% number of wells
state.species0=state.species;%mass of each seven species(H2O,CO2(g),NaCl(s),Na+,Cl,HYD(s),CO2(aq))
state.frac0=state.frac; % mass fractions of four components(water, CO2, Cl, Z) in liquid phase
state.s0=state.s;% saturation of liquid and gas phases (Sl,Sg)
state.ss0=state.ss;% mass fraction of hydrate and halite in the solid phase
state.Tk0=state.Tk; % temperature in Kelvin
state.poro0=state.poro;% porosity (i.e., volume fraction of both fluid phases)
state.perm0=state.perm;% intrinsic permeability
state.mu0=state.mu;% viscosity of liquid and gas phases.
state.rho0=state.rho;% density of liquid and gas phases.
state.pressure0=state.pressure;% liquid and gas pressures.
state.salinity0=state.salinity;% salinity of liquid phase
state.component0=state.component;% masses of four components (water, CO2, Cl, Z) [kg] 
state.m_NaCl0=state.m_NaCl;% molality of aqueous NaCl (we call salinity)
state.poroh0=state.poroh;% volume fraction of hydrate (phi_{HYD})
state.porona0=state.porona;% volume fraction of halite
W=opt.wells;% injection well
for kk=1:nw
    W(kk).pressure0=W(kk).pressure;
end



%% using the chemical reaction module to update the derivative of mass composition with respect to gas pressure and temperture
state.peq=fluid.peq(state);% equilibrium pressure for hydrate formation (interpolated table from temperature and salinity)
STATE=kinetic_reaction_hc(state,G.cells.volumes,fluid,dt,rock);% calculate the chemical reaction to obtain the following derivatives
state.dphiT= STATE.dphiT;% d_porosity/dT
state.dphip= STATE.dphip;%d_porosity/dp_g
state.dspecies= STATE.dspecies;%d_species/dp_g; mass of each species in [kg]
state.dspeciesT= STATE.dspeciesT;%d_species/dT
state.dfrac= STATE.dfrac;% derivitive of mass fraction of each component in water with respect to gas pressure; d_componentfraction/dp_g
state.dfracT= STATE.dfracT;% derivitive of mass fraction of each component in water with respect to temperature; d_componentfraction/dT
state.dsh= STATE.dsh;% derivitive of mass fraction of hydrate and halite in solid phase with respect to gas pressure; 
state.dshT= STATE.dshT;% derivitive of mass fraction of hydrate and halite in solid phase with respect to temperature; 
state.dsp= STATE.dsp;%derivitive of salinity with respect to gas pressure
state.dspT= STATE.dspT;%derivitive of salinity with respect to temperature
state.dcp= STATE.dcp;%derivitive of mass fraction of CO2(aq) in water with respect to gas pressure
state.dcpT= STATE.dcpT;%derivitive of mass fraction of CO2(aq) in water with respect to temperature
state0=state;% save state variable in a temporary space
[d_rho,d_mu]=D_RHOMU_kinetic_h(state);% derivative of density and viscosity over pressures and temperature
[state.h0,state.u0,~,~]=heat_h(state,d_rho);% specific enthalpy and specific internal energy
%% the main loop for reactive transport
%error control 
conlim=1e118;
dplim=state.dplim;% tolerance of error
KNM=20;% maximum iteration number
dp=100000.1;% initial given residual error (this value can be any as long as it is greater than dplim) 
knr=1;% iteration number counter
condL=1;% condition number of Jacobian matrix
realcon=true;% flag for real solution; if realcon=0, then the solution is not real
while (max(abs(dp))>dplim) &&knr<=KNM && condL< conlim &&  realcon


    dp0=dp;% maximum pressure change in the last iteration
%{ 
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
 residual vector R=[-fw;-fc,-fe,-fW] 
 here,fw=governing mass balance equation of water component
      fc=governing mass balance equation of CO2 component
      fHYD=governing mass balance equation of hydrate component
      fe=governing energy balance equation  
      fW=governing mass balance equation of CO2 component in the injection well
      pl=liquid pressure (1st independent transport variable)
      pg=gas pressure (2nd independent transport variable)
      phi=volume fraction of the fluid phase (3rd independent transport variable)
      T=temperature (4th independent transport variable)
      pbh=botttom hole pressure (additional independent transport variable)
 'W' storing all parameters of the injection well  
 'CFL' contains the maximum allowed time step based on CFL condition
 'outflux' contains the mass flux of (water and CO2 components) going out of the boundary
 'flux' stores the Darcy discharge of each face; it can be used to plot the flow field 
 'COMPONENT' stores the masses of four components(water, CO2, Cl, Z) in each cell 
%}

    [L,R,W,CFL,outflux,flux,COMPONENT,~]=LRHS(state,d_rho,d_mu,fluid,G,rock,fht,dt,opt.bc,W);
    condL=condest(L);% condition number of the Jacobian matrix
    state.outflux=outflux;% out mass flux of  water and CO2 components through the open boundary
    state.flux=flux;% Darcy flux of each internel phase
    dx=L\R;% solve the Newton equation to obtain the increase of (p_l,p_g,phi,T) in one (Newton-Raphson) NR iteration;dx=[dp_l,dp_g,dphi,dT]
    nanb=isnan( sum(dx));% check if there is NaN
    msg='Matrix is singular to working precision.';
    msg1=lastwarn;
    if strcmp(msg,msg1) || nanb
         % if there is NaN, then convergency fails
        condL= conlim*2;
        lastwarn('');
    end
    state.CFL=CFL;%CFL condition for updating timestep
    state.pressure(:,1)=state.pressure(:,1)+dx(1:ncc);%p_l=p_l+dp_l
    state.pressure(:,2)=state.pressure(:,2)+dx(1+ncc:2*ncc);%p_g=p_g+dp_g
    state.poro=state.poro+dx(1+2*ncc:3*ncc);%phi=phi+dphi;(volume fraction of both fluid phases)
    state.Tk=state.Tk+dx(1+3*ncc:4*ncc);% T=T+dT
    dpb=dx(1:ncc);% liquid pressure change: dpl
    state.pressure(state.pressure(:,1)<0,1)= (state.pressure(state.pressure(:,1)<0,1)-dpb(state.pressure(:,1)<0,1)).*0.8;% if liquid pressure is lower than zero, pl=pl*0.8
    dpg=dx(ncc+1:2*ncc);% gas pressure change: dpg
    state.pressure(state.pressure(:,2)<0,2)= (state.pressure(state.pressure(:,2)<0,2)-dpg(state.pressure(:,2)<0,1)).*0.8;% if gas pressure is lower than zero, pg=pg*0.8
    state.s=fluid.S(state);% update saturation using the updated pressure
    state.component=COMPONENT;% update the component, which has been explicitly updated inside 'LRHS'

    for k=1:nw
        W(k).pressure=W(k).pressure+dx(4*ncc+k); % update the bottom hole pressure of each well.
    end
    dp=max(abs(dx(1:(4*ncc))));% update the maximum change of liquid pressure, gas pressure, porosity and temperature
     %dp=max(abs(dx));
    realcon=isreal (state.pressure)&&min(state.pressure(:))>0&&min(state.Tk(:))>0&&(~isnan(sum(state.pressure(:))))&&(max(abs(state.pressure(:)))<1e8);% check if there is unreasonable value.
    if realcon&&~nanb
    % perform chemical reaction calculation with the chemical reaction module 'kinetic_reaction_hc'
        state=kinetic_reaction_hc(state,G.cells.volumes,fluid,dt,rock);% update the mass of all species and related derivatives (e.g., d X_l^C/dp) using chemical reaction module
        % porosity is adjusted based on the updated hydrate mass (the 6th column of 'species' contains the mass of HYD species, and the 3rd column of 'rhoref' contains the HYD density)
        state.poro=rock.poro-state.species(:,6)./state.rhoref(3)./G.cells.volumes;
        if isnan(sum(state.species(:,1)))
            state.solve=0;% indicates that the solotion fails.The Newton loop jumps out and restart with a new Newton loop with smaller time step.
            realcon=0;% if we get NaN, then the solution is not reasonable. The Newton loop jumps out and restart with a new Newton loop with smaller time step.
        else
            state.peq=fluid.peq(state); % update the equilibrium pressure for hydrate formation
            [state.rho,state.mu]=rhomu_p_frac_kinetic_h(state);% update the density and viscosity
            [d_rho,d_mu]=D_RHOMU_kinetic_h(state);% update the derivative of density and viscosity
            state.solve=1;% the solution has not failed yet
        end
    end
    knr=knr+1; % count the iteration number
    if dp0>dp && knr>5
        dplim=state.dplim*state.dplim_relax;% loose the convergency control after 5 iterations, provided that dp<dp0
    end

end


%%
if max(abs(state.m_NaCl-state.m_NaCl0))>70e-4
    realcon=0;% the maximum salinity change is limited to 7e-4
end
%% if the previous newton loop fails, we start a new Newton loop with smaller time step
if knr>KNM || condL> conlim ||(~realcon)|| (max(abs(dp))>dplim)||isnan(dp)||nanb
    dplim=state.dplim/100;% reduce the tolerance
    dt=state.meandt*1e-4;% reduce the time step; here 'meandt' is the mean time step in the previous calculations
    dt=min(dt,0.1);% new time step is smaller than 0.1 second
    state=state0;% restart from initial state
    state.species=state.species0; %reinitialize the mass of the species
    [d_rho,d_mu]=D_RHOMU_kinetic_h(state);% deriavtive of density and viscosity
    [state.h0,state.u0,~,~]=heat_h(state,d_rho);% specific enthalphy and specific internal energy
    for kk=1:nw
        W(kk).pressure=W(kk).pressure0;% recover the well pressure
    end
    dp=100; % initial error; and be any value larger than 'dplim'
    %% the following calculation of reactive flow based on Newton method is very similar to the previous section
    while (max(abs(dp))>dplim) &&knr<60
        %obtain the left-hand-side (LHS) matrix and the right-hand-side (RHS) vector.
        [L,R,W,CFL,outflux,flux,COMPONENT,accum]=LRHS(state,d_rho,d_mu,fluid,G,rock,fht,dt,opt.bc,W);
        condL=condest(L);
        state.outflux=outflux;
        state.flux=flux;
        dx=L\R;
        nanb=isnan( sum(dx))||~isreal(sum(dx));
        msg='Matrix is singular to working precision.';
        msg1=lastwarn;
        if strcmp(msg,msg1) || nanb
            condL=conlim*2;
            lastwarn('');
        end
        state.CFL=CFL;
        realcon=isreal (state.pressure);
        if realcon&&~nanb &&(condL< conlim)&&knr<59
            state.pressure(:,1)=state.pressure(:,1)+dx(1:ncc);
            state.pressure(:,2)=state.pressure(:,2)+dx(1+ncc:2*ncc);
            state.poro=state.poro+dx(1+2*ncc:3*ncc);
            state.Tk=state.Tk+dx(1+3*ncc:4*ncc);
            dpb=dx(1:ncc);
            state.pressure(state.pressure(:,1)<0,1)= (state.pressure(state.pressure(:,1)<0,1)-dpb(state.pressure(:,1)<0,1)).*0.8;
            dpg=dx(ncc+1:2*ncc);
            state.pressure(state.pressure(:,2)<0,2)= (state.pressure(state.pressure(:,2)<0,2)-dpg(state.pressure(:,2)<0,1)).*0.8;
            state.s=fluid.S(state);
            state.component=COMPONENT;
            for k=1:nw
                W(k).pressure=W(k).pressure+dx(4*ncc+k);
            end
            dp=(abs(R(1:4*ncc)./accum(:)));
            dp=max(dp(1:2*ncc));
            dp=max(dp,max(abs(dx))./1e2);

            if min(state.pressure(:))<0||min(state.Tk(:))<0||isnan(dp)||~isreal(dx)
                state.solve=0;
                knr=59;
            else
                state.peq=fluid.peq(state);
                if any(isnan(state.peq))||~isreal(state.peq)
                    % pause
                    knr=59;
                end

                state=kinetic_reaction_hc(state,G.cells.volumes,fluid,dt,rock);
                [d_rho,d_mu]=D_RHOMU_kinetic_h(state);
                [state.rho,state.mu]=rhomu_p_frac_kinetic_h(state);
                state.solve=1;
                knr=knr+1;
            end


        else
            aa=randn(1,1)*10;
            state.pressure0(:,2)=state.pressure0(:,2)+aa;
            state.pressure(:,1)=state.pressure0(:,1)+aa;
            state.pressure(:,2)=state.pressure0(:,2)+aa;
            state.s0=fluid.S(state);
            state.Tk=state.Tk0;
            state.h=state.h0;
            state.u=state.u0;
            state.perm=state.perm0;
            state.species=state.species0;
            state.frac=state.frac0;
            state.s=state.s0;
            state.poro=state.poro0;
            state.porona=state.porona0;
            state.poronh=state.poroh0;
            state.ss=state.ss0;
            state.mu=state.mu0;
            state.rho=state.rho0;
            state.pressure=state.pressure0;
            state.salinity=state.salinity0;
            state.component=state.component0;
            state.m_NaCl=state.m_NaCl0;
            state.peq=fluid.peq(state);
            state=kinetic_reaction_hc(state,G.cells.volumes,fluid,dt,rock);
            [d_rho,d_mu]=D_RHOMU_kinetic_h(state);
            [state.h0,state.u0,~,~]=heat_h(state,d_rho);
            state.h=state.h0;
            state.u=state.u0;
            for k=1:nw
                W(k).pressure=W(k).pressure0+abs(randn(1,1)*10);
            end


            [state.rho,state.mu]=rhomu_p_frac_kinetic_h(state);
            dp=dplim*0.1;
            dt=rand(1,1)*state.meandt;
            state.CFL=1e5;
            state.solve=0;

        end






    end
end
%% we have finished the reactive flow calculation, we record the pressure change and saturation change after one time step

dcell=state.pressure-state.pressure0;
dcell=dcell(:);
dp1=zeros(2*ncc+nw,1);
dp1(1:2*ncc)=dcell;
for k=1:nw
    dwell=W(k).pressure(:,1)-W(k).pressure0(:,1);
    dp1(2*ncc+k)=dwell;

end
state.dp(:,1)=dp1;
state.ds=state.s-state.s0;
% the maximum pressure change and saturation change is going to be used to update the next time step 
end









function [L,R,W,CFL,outflux,flux,COMPONENT,accum]=LRHS(state,drho,dmu,fluid,G,rock,fht,dt,bc,W)
%{
 This function serves to calculate the Jacobian matrix and the residual
 vector of the governing flow equations
%} 
% the two ajacent cells of a face; the indicator of the boundary cell is 0. 
[neighborship, ~] = getNeighbourship(G, 'Topological', true);
% flag of the interface
i  = all(neighborship ~= 0, 2);
% the left and right neighbor cells of each interface
inter=neighborship(i,:);
% number of cells
ncc=G.cells.num; 
% saturations state.s=(S_l,S_g) and the derivative of S_l over p_g,
% i.e.,dS=dS_l/dp_g
[state.s,dS]=fluid.S(state);
% relative liquid and gas permeabilities and their derivatives;
% kr=[krl,krg];dkr=[dkrl/dS_l,dkrg/dS_l]
[kr,dkr]=fluid.relperm(state.s);
% get the upwind cell of an internal face; ib= the upwind cell for water
% flow;ic=the upwind cell for the CO2 flow
[ib,ic]=upwind(G,state);
% number of faces; if we have two cubic cells, we will have 11 faces.
nf=G.faces.num;
tib_CFL=fht(i).*kr(ib,1)./state.mu(ib,1).*...
    (state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)+0.5.*(state.rho(inter(:,1),1)+state.rho(inter(:,2),1)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'))...
    ./(G.cells.volumes(inter(:,1)).*rock.poro(inter(:,1))+G.cells.volumes(inter(:,2)).*rock.poro(inter(:,2))).*2;% CFL condition due to liquid phase flow
tic_CFL=fht(i).*kr(ic,2)./state.mu(ic,2).*...
    (state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)+0.5.*(state.rho(inter(:,1),2)+state.rho(inter(:,2),2)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'))...
    ./(G.cells.volumes(inter(:,1)).*rock.poro(inter(:,1))+G.cells.volumes(inter(:,2)).*rock.poro(inter(:,2))).*2;% CFL condition due to gas phase flow
ti_CFL=max([abs(tib_CFL);abs(tic_CFL)]);
CFL=1/ti_CFL;
state.flux=zeros(nf,2);% flux of each phase at each face
tib_CFL=fht(i).*kr(ib,1)./state.mu(ib,1)./G.faces.areas(i).*...
    (state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)+0.5.*(state.rho(inter(:,1),1)+state.rho(inter(:,2),1)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'));
tic_CFL=fht(i).*kr(ic,2)./state.mu(ic,2)./G.faces.areas(i).*...
    (state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)+0.5.*(state.rho(inter(:,1),2)+state.rho(inter(:,2),2)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'));
state.flux(i,1)=tib_CFL;% flux of liquid phase
state.flux(i,2)=tic_CFL;% flux of gas phase
flux=state.flux;% output the flux

%% prepare the advection and diffusion terms
%rho*fht/mu*(dkr/dS) for each face; 'rho' is the density of upwind cell; 'fht' 
%is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'dkr/dS' the derivitive of relative permeability with respect to saturation 
%of the upwind cell; first coulum is for water phase and the second column for CO2 phase
rhokmudkr=zeros(nf,2);
rhokmudkr(i,1)=state.rho(ib,1).*fht(i).*dkr(ib,1)./state.mu(ib,1);
rhokmudkr(i,2)=state.rho(ic,2).*fht(i).*dkr(ic,2)./state.mu(ic,2);

%rho*rho_m*fht/mu*(dkr/dS) for each face; 'rho' is the density of upwind cell; 
%'rho_m' is the mean density of two ajacent cells; 'fht' 
%is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'dkr/dS' the derivitive of relative permeability with respect to saturation 
%of the upwind cell; first coulum is for water phase and the second column for CO2 phase
rhoskmudkr=zeros(nf,2);
rhoskmudkr(i,1)=state.rho(ib,1).*0.5.*(state.rho(inter(:,1),1)+state.rho(inter(:,2),1)).*fht(i).*dkr(ib,1)./state.mu(ib,1);
rhoskmudkr(i,2)=state.rho(ic,2).*0.5.*(state.rho(inter(:,1),2)+state.rho(inter(:,2),2)).*fht(i).*dkr(ic,2)./state.mu(ic,2);

%rho*fht/mu*kr for each face; 'rho' is the density of upwind cell; 'fht' 
%is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'kr' the  relative permeability of the upwind cell;
%first coulum is for water phase and the second column for CO2 phase
rhokmukr=zeros(nf,2);
rhokmukr(i,1)=state.rho(ib,1).*fht(i).*kr(ib,1)./state.mu(ib,1);
rhokmukr(i,2)=state.rho(ic,2).*fht(i).*kr(ic,2)./state.mu(ic,2);

%rho*fht*kr/mu^2*(-dmu/dp) for each face; 'rho' is the density of upwind cell; 'fht' 
%is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'dmu/dp' the derivitive of viscosity with respect to pressure; 'kr' the relative permeability 
%of the upwind cell; first coulum is for water phase and the second column for CO2 phase
rhokkrdmu=zeros(nf,2);
rhokkrdmu(i,1)=state.rho(ib,1).*fht(i).*kr(ib,1)./(state.mu(ib,1)).^2.*(-dmu(ib,1));
rhokkrdmu(i,2)=state.rho(ic,2).*fht(i).*kr(ic,2)./(state.mu(ic,2)).^2.*(-dmu(ic,2));

%rho*rho_m*fht*kr/mu^2*(-dmu/dp) for each face; 'rho' is the density of upwind cell; 
%'rho_m' is the mean density of two ajacent cells; 'fht' 
%is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'dmu/dp' the derivitive of viscosity with respect to pressure; 'kr' the relative permeability 
%of the upwind cell; first coulum is for water phase and the second column for CO2 phase
rhoskkrdmu=zeros(nf,2);
rhoskkrdmu(i,1)=state.rho(ib,1).*0.5.*(state.rho(inter(:,1),1)+state.rho(inter(:,2),1)).*fht(i).*kr(ib,1)./(state.mu(ib,1)).^2.*(-dmu(ib,1));
rhoskkrdmu(i,2)=state.rho(ic,2).*0.5.*(state.rho(inter(:,1),2)+state.rho(inter(:,2),2)).*fht(i).*kr(ic,2)./(state.mu(ic,2)).^2.*(-dmu(ic,2));

%drho*fht*kr/mu for each face; 'drho' is the derivitive of density with respect to pressure of upwind cell; 
% 'fht' is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'kr' the relative permeability of the upwind cell; first coulum is for water phase and the second column for CO2 phase
drhokkrmu=zeros(nf,2);
drhokkrmu(i,1)=drho(ib,1).*fht(i).*kr(ib,1)./state.mu(ib,1);
drhokkrmu(i,2)=drho(ic,2).*fht(i).*kr(ic,2)./state.mu(ic,2);


%rho*fht*kr/mu^2*(-dmu/dT) for each face; 'rho' is the density of upwind cell; 
%'fht' is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'dmu/dT' the derivitive of viscosity with respect to temperature; 'kr' the relative permeability 
%of the upwind cell; first coulum is for water phase and the second column for CO2 phase
rhokkrdmuT=zeros(nf,2);
rhokkrdmuT(i,1)=state.rho(ib,1).*fht(i).*kr(ib,1)./(state.mu(ib,1)).^2.*(-dmu(ib,5));
rhokkrdmuT(i,2)=state.rho(ic,2).*fht(i).*kr(ic,2)./(state.mu(ic,2)).^2.*(-dmu(ic,6));

%rho*rho_m*fht*kr/mu^2*(-dmu/dT) for each face; 'rho' is the density of upwind cell; 
%'rho_m' is the mean density of two ajacent cells; 'fht' 
%is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'dmu/dT' the derivitive of viscosity with respect to temperature; 'kr' the relative permeability 
%of the upwind cell; first coulum is for water phase and the second column for CO2 phase
rhoskkrdmuT=zeros(nf,2);
rhoskkrdmuT(i,1)=state.rho(ib,1).*0.5.*(state.rho(inter(:,1),1)+state.rho(inter(:,2),1)).*fht(i).*kr(ib,1)./(state.mu(ib,1)).^2.*(-dmu(ib,5));
rhoskkrdmuT(i,2)=state.rho(ic,2).*0.5.*(state.rho(inter(:,1),2)+state.rho(inter(:,2),2)).*fht(i).*kr(ic,2)./(state.mu(ic,2)).^2.*(-dmu(ic,6));

%drho_T*rho_m*fht*kr/mu for each face; 'drho_T' is the derivitive of density with respect to temperature of upwind cell; 
%'rho_m' is the mean density of two ajacent cells; 'fht' 
%is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'kr' the relative permeability of the upwind cell; first coulum is for water phase and the second column for CO2 phase
drhokkrmuT=zeros(nf,2);
drhokkrmuT(i,1)=drho(ib,5).*fht(i).*kr(ib,1)./state.mu(ib,1);
drhokkrmuT(i,2)=drho(ic,6).*fht(i).*kr(ic,2)./state.mu(ic,2);

%drho*rho_m*fht*kr/mu for each face; 'drho' is the derivitive of density with respect to pressure of upwind cell; 
%'rho_m' is the mean density of two ajacent cells; 'fht' 
%is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'kr' the relative permeability of the upwind cell; first coulum is for water phase and the second column for CO2 phase
drhoskkrmu=zeros(nf,2);
drhoskkrmu(i,1)=0.5.*(state.rho(inter(:,1),1)+state.rho(inter(:,2),1)).*drho(ib,1).*fht(i).*kr(ib,1)./state.mu(ib,1);
drhoskkrmu(i,2)=0.5.*(state.rho(inter(:,1),2)+state.rho(inter(:,2),2)).*drho(ic,2).*fht(i).*kr(ic,2)./state.mu(ic,2);

%drho_T*rho_m*fht*kr/mu^2*(-dmu/dT) for each face; 'drho_T' is the derivitive of density with respect to temperature of upwind cell; 
%'rho_m' is the mean density of two ajacent cells; 'fht' 
%is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
% 'kr' the relative permeability 
%of the upwind cell; first coulum is for water phase and the second column for CO2 phase
drhoskkrmuT=zeros(nf,2);
drhoskkrmuT(i,1)=0.5.*(state.rho(inter(:,1),1)+state.rho(inter(:,2),1)).*drho(ib,5).*fht(i).*kr(ib,1)./state.mu(ib,1);
drhoskkrmuT(i,2)=0.5.*(state.rho(inter(:,1),2)+state.rho(inter(:,2),2)).*drho(ic,6).*fht(i).*kr(ic,2)./state.mu(ic,2);

%rho*fht*kr/mu^2*(-dmu_C/dp) for each face; 'rho' is the density of upwind cell; 
%'fht' is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'dmu_C/dp' the derivitive of water viscosity with respect to gas pressure and the 
%derivitive of gas viscosity with respect to liquid pressure of the upwind cell; 'kr' the relative permeability 
%of the upwind cell; first coulum is for water phase and the second column for CO2 phase
rhokkrdmuC=zeros(nf,2);
rhokkrdmuC(i,1)=state.rho(ib,1).*fht(i).*kr(ib,1)./(state.mu(ib,1)).^2.*(-dmu(ib,3));
rhokkrdmuC(i,2)=state.rho(ic,2).*fht(i).*kr(ic,2)./(state.mu(ic,2)).^2.*(-dmu(ic,4));

%rho*rho_m*fht*kr/mu^2*(-dmu/dp) for each face; 'rho' is the density of upwind cell; 
%'rho_m' is the mean density of two ajacent cells; 'fht' 
%is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'dmu/dp' the derivitive of water viscosity with respect to gas pressure and the 
%derivitive of gas viscosity with respect to liquid pressure of the upwind cell; 'kr' the relative permeability 
%of the upwind cell; first coulum is for water phase and the second column for CO2 phase
rhoskkrdmuC=zeros(nf,2);
rhoskkrdmuC(i,1)=state.rho(ib,1).*0.5.*(state.rho(inter(:,1),1)+state.rho(inter(:,2),1)).*fht(i).*kr(ib,1)./(state.mu(ib,1)).^2.*(-dmu(ib,3));
rhoskkrdmuC(i,2)=state.rho(ic,2).*0.5.*(state.rho(inter(:,1),2)+state.rho(inter(:,2),2)).*fht(i).*kr(ic,2)./(state.mu(ic,2)).^2.*(-dmu(ic,4));

%drho_C*fht*kr/mu for each face; 'drho_C' is the derivitive of liquid density of upwind cell with respect to gas 
%pressure and the derivitive of gas density of upwind cell with respect to liquid pressure; 
% 'fht' is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'kr' the relative permeability of the upwind cell; first coulum is for water phase and the second column for CO2 phase
drhokkrmuC=zeros(nf,2);
drhokkrmuC(i,1)=drho(ib,3).*fht(i).*kr(ib,1)./state.mu(ib,1);
drhokkrmuC(i,2)=drho(ic,4).*fht(i).*kr(ic,2)./state.mu(ic,2);

%drho_C*rho_m*fht*kr/mu for each face; 'drho_C' is the derivitive of liquid density of upwind cell with respect to gas 
%pressure and the derivitive of gas density of upwind cell with respect to liquid pressure; 
%'rho_m' is the mean density of two ajacent cells; 'fht' 
%is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'kr' the relative permeability of the upwind cell; first coulum is for water phase and the second column for CO2 phase
drhoskkrmuC=zeros(nf,2);
drhoskkrmuC(i,1)=0.5.*(state.rho(inter(:,1),1)+state.rho(inter(:,2),1)).*drho(ib,3).*fht(i).*kr(ib,1)./state.mu(ib,1);
drhoskkrmuC(i,2)=0.5.*(state.rho(inter(:,1),2)+state.rho(inter(:,2),2)).*drho(ic,4).*fht(i).*kr(ic,2)./state.mu(ic,2);

%rho*rho_m*fht*kr/mu for each face; 'rho' is the density of upwind cell; 
%'rho_m' is the mean density of two ajacent cells; 'fht' 
%is the transmissivity of the face; 'mu' is the viscosity of the upwind cell; 
%'kr' the relative permeability of the upwind cell; first coulum is for water phase and the second column for CO2 phase
rhoskmukr=zeros(nf,2);
rhoskmukr(i,1)=0.5.*(state.rho(inter(:,1),1)+state.rho(inter(:,2),1)).*state.rho(ib,1).*fht(i).*kr(ib,1)./state.mu(ib,1);
rhoskmukr(i,2)=0.5.*(state.rho(inter(:,1),2)+state.rho(inter(:,2),2)).*state.rho(ic,2).*fht(i).*kr(ic,2)./state.mu(ic,2);

% calculate the half transmissivity of each face due to mass disperion in water (kb),
% mass dispersion in CO2 (kc), and heat dispersion (kT)
[kb,kc,KT]=DISP_H(state,G,rock);
% mass dispersion transmissivity in water (i.e.,brine)
kb=1./accumarray(G.cells.faces(:,1), 1./kb);
% mass dispersion transmissivity in CO2
kc=1./accumarray(G.cells.faces(:,1), 1./kc);
% heat dispersion transmissivity
KT=1./accumarray(G.cells.faces(:,1), 1./KT);

%dS_m*phi*rho_m*D for each face; 'dS_m' is 1/(S_1+S_2)^2; 
%'phi' is the mean porosity of two ajacent cells; 
%'rho_m' is the mean density of two ajacent cells; 'D' 
%is the DISPERSION transmissivity of the face; first coulum is for water phase and the second column for CO2 phase
dsphirhod=zeros(nf,2);
dsphirhod(i,1)=1./(state.s(inter(:,1),1)+state.s(inter(:,2),1)).^2.*0.5.*(state.rho(inter(:,1),1)+state.rho(inter(:,2),1)).*0.5.*(state.poro(inter(:,1))+state.poro(inter(:,2))).*kb(i);
dsphirhod(i,2)=1./(state.s(inter(:,1),2)+state.s(inter(:,2),2)).^2.*0.5.*(state.rho(inter(:,1),2)+state.rho(inter(:,2),2)).*0.5.*(state.poro(inter(:,1))+state.poro(inter(:,2))).*kc(i);
dsphirhod(isnan(dsphirhod))=0;
dsphirhod(dsphirhod==inf)=0;

%S_m*phi*rho_m*D for each face; 'S_m' is the harmonic mean saturation of two ajacent cells; 
%'phi' is the mean porosity of two ajacent cells; 
%'rho_m' is the mean density of two ajacent cells; 'D' 
%is the DISPERSION transmissivity of the face; first coulum is for water phase and the second column for CO2 phase
sphirhod=zeros(nf,2);
sphirhod(i,1)=2./(1./state.s(inter(:,1),1)+1./state.s(inter(:,2),1)).*0.5.*(state.rho(inter(:,1),1)+state.rho(inter(:,2),1)).*0.5.*(state.poro(inter(:,1))+state.poro(inter(:,2))).*kb(i);
sphirhod(i,2)=2./(1./state.s(inter(:,1),2)+1./state.s(inter(:,2),2)).*0.5.*(state.rho(inter(:,1),2)+state.rho(inter(:,2),2)).*0.5.*(state.poro(inter(:,1))+state.poro(inter(:,2))).*kc(i);
sphirhod(isnan(sphirhod))=0;
sphirhod(sphirhod==inf)=0;%

%S_m*rho_m*D for each face; 'S_m' is the harmonic mean saturation of two ajacent cells; 
%'rho_m' is the mean density of two ajacent cells; 'D' 
%is the DISPERSION transmissivity of the face; first coulum is for water phase and the second column for CO2 phase
srhod(i,1)=2./(1./state.s(inter(:,1),1)+1./state.s(inter(:,2),1)).*0.5.*(state.rho(inter(:,1),1)+state.rho(inter(:,2),1)).*kb(i);
srhod(i,2)=2./(1./state.s(inter(:,1),2)+1./state.s(inter(:,2),2)).*0.5.*(state.rho(inter(:,1),2)+state.rho(inter(:,2),2)).*kc(i);
srhod(isnan(sphirhod))=0;
srhod(sphirhod==inf)=0;%

%S_m*phi*D for each face; 'S_m' is the harmonic mean saturation of two ajacent cells; 
%'phi' is the mean porosity of two ajacent cells; 'D' 
%is the DISPERSION transmissivity of the face; first coulum is for water phase and the second column for CO2 phase
sphid=zeros(nf,2);
sphid(i,1)=2./(1./state.s(inter(:,1),1)+1./state.s(inter(:,2),1)).*0.5.*(state.poro(inter(:,1))+state.poro(inter(:,2))).*kb(i);
sphid(i,2)=2./(1./state.s(inter(:,1),2)+1./state.s(inter(:,2),2)).*0.5.*(state.poro(inter(:,1))+state.poro(inter(:,2))).*kc(i);
sphid(isnan(sphid))=0;
sphid(sphid==inf)=0;%





%% kinetic rate
% equilibrium pressure for given temperature
state.peq=fluid.peq(state);
% flag for HYD formation or dissociation
h_form=state.pressure(:,2)>state.peq;
% reaction area: V*phi^1.5/(2*k)^0.5*r_k0
AS=G.cells.volumes.*state.poro0.^1.5.*(2.*state.permr).^(-.5).*state.kinetic_rate;
% residual molar abundance for chemical reaction; reaction is stopped if
% the molar abundance of the reactant is smaller than the residual value
nr_chemical=state.nr_chemical;
in_res=(state.species(:,1)./0.018./6<nr_chemical|state.species(:,2)./0.044<nr_chemical) & state.pressure(:,2)>state.peq(:);
in_res=in_res|(state.species(:,6)./0.152<nr_chemical& state.pressure(:,2)<state.peq(:));
water_limit=state.species(:,1)./0.018./6<state.species(:,2)./0.044;
AS(in_res)=0;
% effective reaction area
AS1(h_form&water_limit)=AS(h_form&water_limit).*(state.species(h_form&water_limit,1)-nr_chemical(h_form&water_limit).*6.*0.018)./state.rock_mass(h_form&water_limit);
AS1(h_form&~water_limit)=AS(h_form&~water_limit).*(state.species(h_form&~water_limit,2)-nr_chemical(h_form&~water_limit).*0.044)./state.rock_mass(h_form&~water_limit);
AS1(~h_form)=AS(~h_form).*(state.species(~h_form,6)-nr_chemical(~h_form).*0.152)./(state.rock_mass(~h_form)+state.species(~h_form,3)+state.species(~h_form,6));
%activation energy
Ea1=state.Ea1;
% gas constant
RG=8.314;
% kinetic reaction rate
r_k=AS1'.*exp(-Ea1./RG./state.Tk).*(state.pressure(:,2)-state.peq).*0.152;

%% the following calculate dr_k/dp_l, dr_k/dp_g, dr_k/dphi and dr_k/dT
das1_pb=zeros(size(AS1));
das1_pb(h_form&water_limit)=AS(h_form&water_limit)./state.rock_mass(h_form&water_limit).*...
    G.cells.volumes(h_form&water_limit).*state.poro(h_form&water_limit).*(drho(h_form&water_limit,1).*state.s(h_form&water_limit,1)-dS(h_form&water_limit).*state.rho(h_form&water_limit,1));

das1_pb(h_form&~water_limit)=AS(h_form&~water_limit)./state.rock_mass(h_form&~water_limit).*...
    G.cells.volumes(h_form&~water_limit).*state.poro(h_form&~water_limit).*(dS(h_form&~water_limit).*state.rho(h_form&~water_limit,2));


das1_pg=zeros(size(AS1));
das1_pg(h_form&water_limit)=AS(h_form&water_limit)./state.rock_mass(h_form&water_limit).*...
    G.cells.volumes(h_form&water_limit).*state.poro(h_form&water_limit).*(drho(h_form&water_limit,3).*state.s(h_form&water_limit,1)+dS(h_form&water_limit).*state.rho(h_form&water_limit,1));

das1_pg(h_form&~water_limit)=AS(h_form&~water_limit)./state.rock_mass(h_form&~water_limit).*...
    G.cells.volumes(h_form&~water_limit).*state.poro(h_form&~water_limit).*(drho(h_form&~water_limit,2).*state.s(h_form&~water_limit,2)-dS(h_form&~water_limit).*state.rho(h_form&~water_limit,2));



das1_phi=zeros(size(AS1));
das1_phi(h_form&water_limit)=AS(h_form&water_limit)./state.rock_mass(h_form&water_limit).*...
    G.cells.volumes(h_form&water_limit).*(state.rho(h_form&water_limit,1).*state.s(h_form&water_limit,1));

das1_phi(~h_form)=AS(~h_form).*...
    (state.rhoref(3).*G.cells.volumes(~h_form).*(-1).*(state.rock_mass(~h_form)+state.species(~h_form,3)+state.species(~h_form,6))-...
    state.rhoref(3).*G.cells.volumes(~h_form).*(-1).*(state.species(~h_form,6)))./...
    (state.rock_mass(~h_form)+state.species(~h_form,3)+state.species(~h_form,6)).^2;

dr_k_pb=das1_pb'.*exp(-Ea1./RG./state.Tk).*(state.pressure(:,2)-state.peq).*0.152;
dr_k_pg=das1_pg'.*exp(-Ea1./RG./state.Tk).*(state.pressure(:,2)-state.peq).*0.152+...
    AS1'.*exp(-Ea1./RG./state.Tk).*0.152;
dr_k_phi=das1_phi'.*exp(-Ea1./RG./state.Tk).*(state.pressure(:,2)-state.peq).*0.152;

dT_peq=0.01;
state.Tk=state.Tk+dT_peq;
peq_n=fluid.peq(state);
state.Tk=state.Tk-dT_peq;
dpeq_T=(peq_n-state.peq)./dT_peq;

dr_k_T=AS1'.*exp(-Ea1./RG./state.Tk).*Ea1./RG./state.Tk.^2.*(state.pressure(:,2)-state.peq).*0.152+...
    AS1'.*exp(-Ea1./RG./state.Tk).*(-dpeq_T).*0.152;




%%
% state.peq=fluid.peq(state);
% h_form=state.pressure(:,2)>state.peq;
% AS=G.cells.volumes.*state.poro0.^1.5.*(2.*state.permr).^(-.5).*state.kinetic_rate;
% nr_chemical=state.nr_chemical;
% in_res=(state.species(:,1)./0.018./6<nr_chemical|state.species(:,2)./0.044<nr_chemical) & state.pressure(:,2)>state.peq(:);
% in_res=in_res|(state.species(:,6)./0.152<nr_chemical& state.pressure(:,2)<state.peq(:));
% AS(in_res)=0;
%
%
%
%
%
% water_limit=state.species(:,1)./0.018./6<state.species(:,2)./0.044;
% AS1=zeros(size(AS));
% AS1(h_form&water_limit)=AS(h_form&water_limit).*(state.species(h_form&water_limit,1)-nr_chemical(h_form&water_limit).*6.*0.018)./state.rock_mass(h_form&water_limit);
% AS1(h_form&~water_limit)=AS(h_form&~water_limit).*(state.species(h_form&~water_limit,2)-nr_chemical(h_form&~water_limit).*0.044)./state.rock_mass(h_form&~water_limit);
% AS1(~h_form)=AS(~h_form).*(state.species(~h_form,6)-nr_chemical(~h_form).*0.152)./(state.rock_mass(~h_form));
% Ea1=state.Ea1;
% RG=8.314;
% r_k=AS1.*exp(-Ea1./RG./state.Tk).*(state.pressure(:,2)-state.peq).*0.152;
%
%
%
% das1_pb=zeros(size(AS1));
% das1_pb(h_form&water_limit)=AS(h_form&water_limit)./state.rock_mass(h_form&water_limit).*...
%     G.cells.volumes(h_form&water_limit).*state.poro(h_form&water_limit).*(drho(h_form&water_limit,1).*state.s(h_form&water_limit,1)-dS(h_form&water_limit).*state.rho(h_form&water_limit,1));
%
% das1_pb(h_form&~water_limit)=AS(h_form&~water_limit)./state.rock_mass(h_form&~water_limit).*...
%     G.cells.volumes(h_form&~water_limit).*state.poro(h_form&~water_limit).*(dS(h_form&~water_limit).*state.rho(h_form&~water_limit,2));
%
%
% das1_pg=zeros(size(AS1));
% das1_pg(h_form&water_limit)=AS(h_form&water_limit)./state.rock_mass(h_form&water_limit).*...
%     G.cells.volumes(h_form&water_limit).*state.poro(h_form&water_limit).*(drho(h_form&water_limit,3).*state.s(h_form&water_limit,1)+dS(h_form&water_limit).*state.rho(h_form&water_limit,1));
%
% das1_pg(h_form&~water_limit)=AS(h_form&~water_limit)./state.rock_mass(h_form&~water_limit).*...
%     G.cells.volumes(h_form&~water_limit).*state.poro(h_form&~water_limit).*(drho(h_form&~water_limit,2).*state.s(h_form&~water_limit,2)-dS(h_form&~water_limit).*state.rho(h_form&~water_limit,2));
%
%
%
% das1_phi=zeros(size(AS1));
% das1_phi(h_form&water_limit)=AS(h_form&water_limit)./state.rock_mass(h_form&water_limit).*...
%     G.cells.volumes(h_form&water_limit).*(state.rho(h_form&water_limit,1).*state.s(h_form&water_limit,1));
% das1_phi(h_form&~water_limit)=AS(h_form&~water_limit)./state.rock_mass(h_form&~water_limit).*...
%     G.cells.volumes(h_form&~water_limit).*(state.rho(h_form&~water_limit,2).*state.s(h_form&~water_limit,2));
% das1_phi(~h_form)=AS(~h_form).*...
%     (-state.rhoref(3).*G.cells.volumes(~h_form))./...
%     (state.rock_mass(~h_form));
%
% % das1_pb=zeros(size(AS1));
% % das1_pg=zeros(size(AS1));
% % das1_phi=zeros(size(AS1));
% dr_k_pb=0.*das1_pb.*exp(-Ea1./RG./state.Tk).*(state.pressure(:,2)-state.peq).*0.152;
% dr_k_pg=das1_pg.*exp(-Ea1./RG./state.Tk).*(state.pressure(:,2)-state.peq).*0.152+...
%     AS1.*exp(-Ea1./RG./state.Tk).*0.152;
%
%
%
%
% dr_k_pg= AS1.*exp(-Ea1./RG./state.Tk).*0.152;
% dr_k_phi=das1_phi.*exp(-Ea1./RG./state.Tk).*(state.pressure(:,2)-state.peq).*0.152;
%
% dT_peq=-0.01;
% state.Tk=state.Tk+dT_peq;
% peq_n=fluid.peq(state);
% state.Tk=state.Tk-dT_peq;
% dpeq_T=(peq_n-state.peq)./dT_peq;
%
% dr_k_T=AS1.*exp(-Ea1./RG./state.Tk).*Ea1./RG./state.Tk.^2.*(state.pressure(:,2)-state.peq).*0.152+...
%     AS1.*exp(-Ea1./RG./state.Tk).*(-dpeq_T).*0.152;
%
% dr_k_phi=zeros(size(dr_k_phi));
%
%
%
% dr_k_pg=1.*state.dspecies(:,6)./dt;
% dr_k_T=1.*state.dspeciesT(:,6)./dt;
%% m11: df_w/dp_l

%{
 jacobian matrix L=[{dfw/dpl},dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}

% derivitive of storage term for the water component  with respect to
% liquid pressure change 
d1e=G.cells.volumes.*state.poro./dt.*(state.s(:,1).*drho(:,1).*1.*state.frac(:,1)+state.rho(:,1).*(-dS).*state.frac(:,1));
% derivative of the mass change due to  the hydration reaction with respect
% to liquid pressure
d1e=d1e+1.*18.*6./152.*dr_k_pb;

m11=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));
% add the derivative terms due to advection diffusion 
m11=m11...
    +sparse(inter(:,1),ib,state.frac(ib,1).* rhokkrdmu(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* rhokkrdmu(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,1).* rhoskkrdmu(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* rhoskkrdmu(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.frac(ib,1).* drhokkrmu(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* drhokkrmu(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,1).* drhoskkrmu(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* drhoskkrmu(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),state.frac(ib,1).* rhokmukr(i,1).*(-1),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1),state.frac(ib,1).* rhokmukr(i,1).*(1),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.frac(ib,1).* rhokmukr(i,1).*(-1),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.frac(ib,1).* rhokmukr(i,1).*(1),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),state.frac(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,2),1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1),state.frac(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,1),1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.frac(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,1),1).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.frac(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,2),1).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.frac(ib,1).* rhokmudkr(i,1).*(-dS(ib)).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* rhokmudkr(i,1).*(-dS(ib)).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,1).* rhoskmudkr(i,1).*(-dS(ib)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* rhoskmudkr(i,1).*(-dS(ib)).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), sphid(i,1).*0.5.*drho(inter(:,2),1).*(-state.frac(inter(:,2),1)+state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), sphid(i,1).*0.5.*drho(inter(:,1),1).*(-state.frac(inter(:,2),1)+state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), sphid(i,1).*0.5.*drho(inter(:,2),1).*(state.frac(inter(:,2),1)-state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),sphid(i,1).*0.5.*drho(inter(:,1),1).*(state.frac(inter(:,2),1)-state.frac(inter(:,1),1)),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), 2.*state.s(inter(:,1),1).^2.*dsphirhod(i,1).*(-dS(inter(:,2))).*(-state.frac(inter(:,2),1)+state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), 2.*state.s(inter(:,2),1).^2.*dsphirhod(i,1).*(-dS(inter(:,1))).*(-state.frac(inter(:,2),1)+state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), 2.*state.s(inter(:,1),1).^2.*dsphirhod(i,1).*(-dS(inter(:,2))).*(-state.frac(inter(:,1),1)+state.frac(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1), 2.*state.s(inter(:,2),1).^2.*dsphirhod(i,1).*(-dS(inter(:,1))).*(-state.frac(inter(:,1),1)+state.frac(inter(:,2),1)),(ncc),(ncc));




%% m12: dfw/dp_g
%{
 jacobian matrix L=[dfw/dpl,{dfw/dpg},dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}


%derivitive of storage term for the water component with respect to gas
%pressure change 
d1e=G.cells.volumes.*state.poro./dt.*(drho(:,3).*state.s(:,1).*state.frac(:,1)+state.rho(:,1).*dS.*state.frac(:,1)+state.s(:,1).*state.rho(:,1).*state.dfrac(:,1));
% derivative of the mass change due to  the hydration reaction with respect
% to gas pressure
d1e=d1e+1.*18.*6./152.*dr_k_pg;


% add the derivative terms due to advection diffusion 
m12=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));
m12=m12...
    +sparse(inter(:,1),ib,state.frac(ib,1).* rhokkrdmuC(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...%NEW
    +sparse(inter(:,2),ib,state.frac(ib,1).* rhokkrdmuC(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,1).* rhoskkrdmuC(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* rhoskkrdmuC(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.frac(ib,1).* drhokkrmuC(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...%NEW
    +sparse(inter(:,2),ib,state.frac(ib,1).* drhokkrmuC(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,1).* drhoskkrmuC(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* drhoskkrmuC(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),state.frac(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,2),3).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...%NEW
    +sparse(inter(:,1),inter(:,1),state.frac(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,1),3).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.frac(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,1),3).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.frac(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,2),3).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.dfrac(ib,1).* rhokmukr(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.dfrac(ib,1).* rhokmukr(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.dfrac(ib,1).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.dfrac(ib,1).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.frac(ib,1).* rhokmudkr(i,1).*(dS(ib)).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* rhokmudkr(i,1).*(dS(ib)).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,1).* rhoskmudkr(i,1).*(dS(ib)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* rhoskmudkr(i,1).*(dS(ib)).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), sphid(i,1).*0.5.*drho(inter(:,2),3).*(-state.frac(inter(:,2),1)+state.frac(inter(:,1),1)),(ncc),(ncc))...%NEW
    +sparse(inter(:,1),inter(:,1), sphid(i,1).*0.5.*drho(inter(:,1),3).*(-state.frac(inter(:,2),1)+state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), sphid(i,1).*0.5.*drho(inter(:,2),3).*(state.frac(inter(:,2),1)-state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),sphid(i,1).*0.5.*drho(inter(:,1),3).*(state.frac(inter(:,2),1)-state.frac(inter(:,1),1)),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), 2.*state.s(inter(:,1),1).^2.*dsphirhod(i,1).*(-dS(inter(:,2))).*(state.frac(inter(:,2),1)-state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), 2.*state.s(inter(:,2),1).^2.*dsphirhod(i,1).*(-dS(inter(:,1))).*(state.frac(inter(:,2),1)-state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), 2.*state.s(inter(:,1),1).^2.*dsphirhod(i,1).*(-dS(inter(:,2))).*(state.frac(inter(:,1),1)-state.frac(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1), 2.*state.s(inter(:,2),1).^2.*dsphirhod(i,1).*(-dS(inter(:,1))).*(state.frac(inter(:,1),1)-state.frac(inter(:,2),1)),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), -sphirhod(i,1).*(state.dfrac(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), sphirhod(i,1).*(state.dfrac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), sphirhod(i,1).*(state.dfrac(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),-sphirhod(i,1).*(state.dfrac(inter(:,1),1)),(ncc),(ncc));


%% m13: dfw/dphi
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,{dfw/dphi},dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}


%derivitive of storage term for the water component with respect to gas
%pressure change 
d1e=G.cells.volumes./dt.*(state.s(:,1).*state.rho(:,1).*state.frac(:,1));
% derivative of the mass change due to  the hydration reaction with respect
% to porosity
d1e=d1e+1.*18.*6./152.*dr_k_phi;
m13=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));
% add the derivative terms due to advection diffusion 
m13=m13...
    +sparse(inter(:,1),inter(:,2), srhod(i,1).*0.5.*(-state.frac(inter(:,2),1)+state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), srhod(i,1).*0.5.*(-state.frac(inter(:,2),1)+state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), srhod(i,1).*0.5.*(state.frac(inter(:,2),1)-state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),srhod(i,1).*0.5.*(state.frac(inter(:,2),1)-state.frac(inter(:,1),1)),(ncc),(ncc));


%% m14: dfw/dT
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,{dfw/dT},dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}


%derivitive of storage term for the water component with respect to temperature change 

d1e=G.cells.volumes.*state.poro./dt.*(state.s(:,1).*drho(:,5).*state.frac(:,1)+state.s(:,1).*state.rho(:,1).*state.dfracT(:,1));
% derivative of the mass change due to  the hydration reaction with respect
% to temperature
d1e=d1e+1.*18.*6./152.*dr_k_T;
m14=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));
% add the derivative terms due to advection diffusion
m14=m14...
    +sparse(inter(:,1),ib,state.frac(ib,1).* rhokkrdmuT(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* rhokkrdmuT(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,1).* rhoskkrdmuT(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* rhoskkrdmuT(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.dfracT(ib,1).* rhokmukr(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.dfracT(ib,1).* rhokmukr(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.dfracT(ib,1).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.dfracT(ib,1).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.frac(ib,1).* drhokkrmuT(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* drhokkrmuT(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,1).* drhoskkrmuT(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,1).* drhoskkrmuT(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),state.frac(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,2),5).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1),state.frac(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,1),5).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.frac(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,1),5).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.frac(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,2),5).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), sphid(i,1).*0.5.*drho(inter(:,2),5).*(-state.frac(inter(:,2),1)+state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), sphid(i,1).*0.5.*drho(inter(:,1),5).*(-state.frac(inter(:,2),1)+state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), sphid(i,1).*0.5.*drho(inter(:,2),5).*(state.frac(inter(:,2),1)-state.frac(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),sphid(i,1).*0.5.*drho(inter(:,1),5).*(state.frac(inter(:,2),1)-state.frac(inter(:,1),1)),(ncc),(ncc));



%% m21: dfc/dpl
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  {dfc/dpl},dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}


%derivitive of storage term for the CO2 component with respect to liquid pressure change 
d1e=G.cells.volumes.*state.poro./dt.*(state.s(:,1).*drho(:,1).*1.*state.frac(:,2)+state.rho(:,1).*(-dS).*state.frac(:,2)-state.rho(:,2).*(-dS));
% derivative of the mass change due to  the hydration reaction with respect
% to liquid pressure
d1e=d1e+1.*44./152.*dr_k_pb;
% add the derivative terms due to advection diffusion
m21=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));
m21=m21...
    +sparse(inter(:,1),ib,state.frac(ib,2).* rhokkrdmu(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* rhokkrdmu(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,2).* rhoskkrdmu(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* rhoskkrdmu(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.frac(ib,2).* drhokkrmu(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* drhokkrmu(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,2).* drhoskkrmu(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* drhoskkrmu(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),state.frac(ib,2).* rhokmukr(i,1).*(-1),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1),state.frac(ib,2).* rhokmukr(i,1).*(1),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.frac(ib,2).* rhokmukr(i,1).*(-1),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.frac(ib,2).* rhokmukr(i,1).*(1),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),state.frac(ib,2).* rhokmukr(i,1).*0.5.*drho(inter(:,2),1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1),state.frac(ib,2).* rhokmukr(i,1).*0.5.*drho(inter(:,1),1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.frac(ib,2).* rhokmukr(i,1).*0.5.*drho(inter(:,1),1).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.frac(ib,2).* rhokmukr(i,1).*0.5.*drho(inter(:,2),1).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.frac(ib,2).* rhokmudkr(i,1).*(-dS(ib)).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* rhokmudkr(i,1).*(-dS(ib)).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,2).* rhoskmudkr(i,1).*(-dS(ib)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* rhoskmudkr(i,1).*(-dS(ib)).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ic, rhokmudkr(i,2).*(-dS(ic)).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic,rhokmudkr(i,2).*(-dS(ic)).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic, rhoskmudkr(i,2).*(-dS(ic)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic, rhoskmudkr(i,2).*(-dS(ic)).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), sphid(i,1).*0.5.*drho(inter(:,2),1).*(-state.frac(inter(:,2),2)+state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), sphid(i,1).*0.5.*drho(inter(:,1),1).*(-state.frac(inter(:,2),2)+state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), sphid(i,1).*0.5.*drho(inter(:,2),1).*(state.frac(inter(:,2),2)-state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),sphid(i,1).*0.5.*drho(inter(:,1),1).*(state.frac(inter(:,2),2)-state.frac(inter(:,1),2)),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), 2.*state.s(inter(:,1),1).^2.*dsphirhod(i,1).*(-dS(inter(:,2))).*(-state.frac(inter(:,2),2)+state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), 2.*state.s(inter(:,2),1).^2.*dsphirhod(i,1).*(-dS(inter(:,1))).*(-state.frac(inter(:,2),2)+state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), 2.*state.s(inter(:,1),1).^2.*dsphirhod(i,1).*(-dS(inter(:,2))).*(-state.frac(inter(:,1),2)+state.frac(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1), 2.*state.s(inter(:,2),1).^2.*dsphirhod(i,1).*(-dS(inter(:,1))).*(-state.frac(inter(:,1),2)+state.frac(inter(:,2),2)),(ncc),(ncc));



%% m22: dfc/dpg
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,{dfc/dpg},dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}


%derivitive of storage term for the CO2 component with respect to gas pressure change 
d1e=G.cells.volumes.*state.poro./dt.*(drho(:,3).*state.s(:,1).*state.frac(:,2)+state.rho(:,1).*dS.*state.frac(:,2)+state.s(:,1).*state.rho(:,1).*state.dfrac(:,2)...
    -state.rho(:,2).*dS+state.s(:,2).*drho(:,2).*1);
% derivative of the mass change due to  the hydration reaction with respect
% to gas pressure
d1e=d1e+1.*44./152.*dr_k_pg;
m22=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));
% add the derivative terms due to advection diffusion
m22=m22...
    +sparse(inter(:,1),ic,rhokkrdmu(i,2).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic, rhokkrdmu(i,2).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic, rhoskkrdmu(i,2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic, rhoskkrdmu(i,2).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ic, drhokkrmu(i,2).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic, drhokkrmu(i,2).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic,drhoskkrmu(i,2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic, drhoskkrmu(i,2).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),rhokmukr(i,2).*(-1),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), rhokmukr(i,2).*(1),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1), rhokmukr(i,2).*(-1),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), rhokmukr(i,2).*(1),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), rhokmukr(i,2).*0.5.*drho(inter(:,2),2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), rhokmukr(i,2).*0.5.*drho(inter(:,1),2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1), rhokmukr(i,2).*0.5.*drho(inter(:,1),2).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), rhokmukr(i,2).*0.5.*drho(inter(:,2),2).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ic, rhokmudkr(i,2).*(dS(ic)).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic, rhokmudkr(i,2).*(dS(ic)).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic, rhoskmudkr(i,2).*(dS(ic)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic, rhoskmudkr(i,2).*(dS(ic)).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.frac(ib,2).* rhokkrdmuC(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...%NEW
    +sparse(inter(:,2),ib,state.frac(ib,2).* rhokkrdmuC(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,2).* rhoskkrdmuC(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* rhoskkrdmuC(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.frac(ib,2).* drhokkrmuC(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...%NEW
    +sparse(inter(:,2),ib,state.frac(ib,2).* drhokkrmuC(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,2).* drhoskkrmuC(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* drhoskkrmuC(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),state.frac(ib,2).* rhokmukr(i,1).*0.5.*drho(inter(:,2),3).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...%NEW
    +sparse(inter(:,1),inter(:,1),state.frac(ib,2).* rhokmukr(i,1).*0.5.*drho(inter(:,1),3).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.frac(ib,2).* rhokmukr(i,1).*0.5.*drho(inter(:,1),3).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.frac(ib,2).* rhokmukr(i,1).*0.5.*drho(inter(:,2),3).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.dfrac(ib,2).* rhokmukr(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.dfrac(ib,2).* rhokmukr(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.dfrac(ib,2).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.dfrac(ib,2).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.frac(ib,2).* rhokmudkr(i,1).*(dS(ib)).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* rhokmudkr(i,1).*(dS(ib)).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,2).* rhoskmudkr(i,1).*(dS(ib)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* rhoskmudkr(i,1).*(dS(ib)).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), sphid(i,1).*0.5.*drho(inter(:,2),3).*(-state.frac(inter(:,2),2)+state.frac(inter(:,1),2)),(ncc),(ncc))...%NEW
    +sparse(inter(:,1),inter(:,1), sphid(i,1).*0.5.*drho(inter(:,1),3).*(-state.frac(inter(:,2),2)+state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), sphid(i,1).*0.5.*drho(inter(:,2),3).*(state.frac(inter(:,2),2)-state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),sphid(i,1).*0.5.*drho(inter(:,1),3).*(state.frac(inter(:,2),2)-state.frac(inter(:,1),2)),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), 2.*state.s(inter(:,1),1).^2.*dsphirhod(i,1).*(-dS(inter(:,2))).*(state.frac(inter(:,2),2)-state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), 2.*state.s(inter(:,2),1).^2.*dsphirhod(i,1).*(-dS(inter(:,1))).*(state.frac(inter(:,2),2)-state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), 2.*state.s(inter(:,1),1).^2.*dsphirhod(i,1).*(-dS(inter(:,2))).*(state.frac(inter(:,1),2)-state.frac(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1), 2.*state.s(inter(:,2),1).^2.*dsphirhod(i,1).*(-dS(inter(:,1))).*(state.frac(inter(:,1),2)-state.frac(inter(:,2),2)),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), -sphirhod(i,1).*(state.dfrac(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), sphirhod(i,1).*(state.dfrac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), sphirhod(i,1).*(state.dfrac(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),-sphirhod(i,1).*(state.dfrac(inter(:,1),2)),(ncc),(ncc));

%% m23: dfc/dphi
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,{dfc/dphi},dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}


%derivitive of storage term for the CO2 component with respect to porosity change 
d1e=G.cells.volumes./dt.*(state.s(:,1).*state.rho(:,1).*state.frac(:,2)+state.s(:,2).*state.rho(:,2));
% derivative of the mass change due to  the hydration reaction with respect
% to porosity
d1e=d1e+1.*44./152.*dr_k_phi;
m23=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));
% add the derivative terms due to  diffusion
m23=m23...
    ...
    ...
    +sparse(inter(:,1),inter(:,2), srhod(i,1).*0.5.*(-state.frac(inter(:,2),2)+state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), srhod(i,1).*0.5.*(-state.frac(inter(:,2),2)+state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), srhod(i,1).*0.5.*(state.frac(inter(:,2),2)-state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),srhod(i,1).*0.5.*(state.frac(inter(:,2),2)-state.frac(inter(:,1),2)),(ncc),(ncc));


%% m24: dfc/dT
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,{dfc/dT},dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}


%derivitive of storage term for the CO2 component with respect to porosity change 
d1e=G.cells.volumes.*state.poro./dt.*(state.s(:,1).*drho(:,5).*state.frac(:,2)+state.s(:,2).*drho(:,6)+state.s(:,1).*state.rho(:,1).*state.dfracT(:,2));
% derivative of the mass change due to  the hydration reaction with respect
% to temperature
d1e=d1e+1.*44./152.*dr_k_T;
m24=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));
% add the derivative terms due to advection dispersion
m24=m24...
    +sparse(inter(:,1),ib,state.frac(ib,2).* rhokkrdmuT(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* rhokkrdmuT(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,2).* rhoskkrdmuT(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* rhoskkrdmuT(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.frac(ib,2).* drhokkrmuT(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* drhokkrmuT(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.frac(ib,2).* drhoskkrmuT(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.frac(ib,2).* drhoskkrmuT(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),state.frac(ib,2).* rhokmukr(i,1).*0.5.*drho(inter(:,2),5).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1),state.frac(ib,2).* rhokmukr(i,1).*0.5.*drho(inter(:,1),5).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.frac(ib,2).* rhokmukr(i,1).*0.5.*drho(inter(:,1),5).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.frac(ib,2).* rhokmukr(i,1).*0.5.*drho(inter(:,2),5).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), sphid(i,1).*0.5.*drho(inter(:,2),5).*(-state.frac(inter(:,2),2)+state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), sphid(i,1).*0.5.*drho(inter(:,1),5).*(-state.frac(inter(:,2),2)+state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), sphid(i,1).*0.5.*drho(inter(:,2),5).*(state.frac(inter(:,2),2)-state.frac(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),sphid(i,1).*0.5.*drho(inter(:,1),5).*(state.frac(inter(:,2),2)-state.frac(inter(:,1),2)),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ic,rhokkrdmuT(i,2).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic, rhokkrdmuT(i,2).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic, rhoskkrdmuT(i,2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic, rhoskkrdmuT(i,2).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ic, drhokkrmuT(i,2).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic, drhokkrmuT(i,2).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic,drhoskkrmuT(i,2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic, drhoskkrmuT(i,2).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    ...
    +sparse(inter(:,1),inter(:,2), rhokmukr(i,2).*0.5.*drho(inter(:,2),6).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), rhokmukr(i,2).*0.5.*drho(inter(:,1),6).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1), rhokmukr(i,2).*0.5.*drho(inter(:,1),6).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), rhokmukr(i,2).*0.5.*drho(inter(:,2),6).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc));






%% m31: dfHYD/dpl
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  {dfHYD/dpl},dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}


% derivative of the mass change due to  the hydration reaction with respect
% to liquid pressure
d1e=-dr_k_pb;
m31=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));

%% m32: dfHYD/dpg
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,{dfHYD/dpg},dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}


% derivative of the mass change due to  the hydration reaction with respect
% to gas pressure
d1e=-dr_k_pg;
m32=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));



%% m33: dfHYD/dphi
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,{dfHYD/dphi},dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}


%derivitive of storage term for the CO2 component with respect to porosity pressure change 
d1e=-G.cells.volumes.*state.rhoref(3)./dt-dr_k_phi;
m33=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));



%% m34: dfHYD/dT
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,{dfHYD/dT},dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}
% derivative of the mass change due to  the hydration reaction with respect
% to temperature
d1e=-dr_k_T;
m34=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));



%% m41: dfe/dpl
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  {dfe/dpl},dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}
% derivative of the internal energy due to  the hydration reaction with respect
% to liquid pressure 
[state.h,state.u,state.dh,du]=heat_h(state,drho);

heatrate=0;
if heatrate==1
    heatc=1;
else
    heatc=0;
end
% the storage term
d1e=G.cells.volumes.*state.poro./dt.*(state.s(:,1).*drho(:,1).*1.*state.u(:,1)+state.rho(:,1).*(-dS).*state.u(:,1)+state.rho(:,1).*state.s(:,1).*du(:,1)...
    +state.s(:,2).*drho(:,4).*state.u(:,2)+state.rho(:,2).*(dS).*state.u(:,2));
% the chemical reaction terms (not considered, heatc=0)
d1e=d1e-heatc.*dr_k_pb./0.152.*(state.u(:,3).*0.152-state.u(:,2).*0.044-state.u(:,1).*6.*0.018);

m41=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));
% the advection dispersion term
m41=m41...
    +sparse(inter(:,1),ib,state.h(ib,1).* rhokkrdmu(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* rhokkrdmu(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.h(ib,1).* rhoskkrdmu(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* rhoskkrdmu(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.h(ib,1).* drhokkrmu(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* drhokkrmu(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.h(ib,1).* drhoskkrmu(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* drhoskkrmu(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),state.h(ib,1).* rhokmukr(i,1).*(-1),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1),state.h(ib,1).* rhokmukr(i,1).*(1),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.h(ib,1).* rhokmukr(i,1).*(-1),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.h(ib,1).* rhokmukr(i,1).*(1),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),state.h(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,2),1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1),state.h(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,1),1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.h(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,1),1).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.h(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,2),1).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.h(ib,1).* rhokmudkr(i,1).*(-dS(ib)).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* rhokmudkr(i,1).*(-dS(ib)).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.h(ib,1).* rhoskmudkr(i,1).*(-dS(ib)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* rhoskmudkr(i,1).*(-dS(ib)).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.dh(ib,1).* rhokmukr(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.dh(ib,1).* rhokmukr(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.dh(ib,1).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.dh(ib,1).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ic, state.h(ic,2).*rhokmudkr(i,2).*(-dS(ic)).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic,state.h(ic,2).*rhokmudkr(i,2).*(-dS(ic)).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic, state.h(ic,2).*rhoskmudkr(i,2).*(-dS(ic)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic, state.h(ic,2).*rhoskmudkr(i,2).*(-dS(ic)).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc));



%% m42: dfe/dpg
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,{dfe/dpg},dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}
% derivative of the internal energy due to  the hydration reaction with respect
% to gas pressure 
d1e=G.cells.volumes.*state.poro./dt.*(state.s(:,2).*drho(:,2).*1.*state.u(:,2)+state.rho(:,2).*(-dS).*state.u(:,2)+state.rho(:,2).*state.s(:,2).*du(:,2)...
    +state.s(:,1).*drho(:,3).*state.u(:,1)+state.rho(:,1).*(dS).*state.u(:,1));
% the chemical reaction terms (not considered, heatc=0)
d1e=d1e-heatc.*dr_k_pg./0.152.*(state.u(:,3).*0.152-state.u(:,2).*0.044-state.u(:,1).*6.*0.018);
m42=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));
% the advection dispersion term
m42=m42...
    +sparse(inter(:,1),ib,state.h(ib,1).* rhokkrdmuC(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...%NEW
    +sparse(inter(:,2),ib,state.h(ib,1).* rhokkrdmuC(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.h(ib,1).* rhoskkrdmuC(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* rhoskkrdmuC(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.h(ib,1).* drhokkrmuC(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...%NEW
    +sparse(inter(:,2),ib,state.h(ib,1).* drhokkrmuC(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.h(ib,1).* drhoskkrmuC(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* drhoskkrmuC(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),state.h(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,2),3).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...%NEW
    +sparse(inter(:,1),inter(:,1),state.h(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,1),3).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.h(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,1),3).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.h(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,2),3).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.h(ib,1).* rhokmudkr(i,1).*(dS(ib)).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* rhokmudkr(i,1).*(dS(ib)).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.h(ib,1).* rhoskmudkr(i,1).*(dS(ib)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* rhoskmudkr(i,1).*(dS(ib)).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ic,state.h(ic,2).*rhokkrdmu(i,2).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic,state.h(ic,2).* rhokkrdmu(i,2).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic,state.h(ic,2).* rhoskkrdmu(i,2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic, state.h(ic,2).*rhoskkrdmu(i,2).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ic, state.h(ic,2).*drhokkrmu(i,2).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic, state.h(ic,2).*drhokkrmu(i,2).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic,state.h(ic,2).*drhoskkrmu(i,2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic, state.h(ic,2).*drhoskkrmu(i,2).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),state.h(ic,2).*rhokmukr(i,2).*(-1),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1),state.h(ic,2).* rhokmukr(i,2).*(1),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.h(ic,2).* rhokmukr(i,2).*(-1),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.h(ic,2).*rhokmukr(i,2).*(1),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), state.h(ic,2).*rhokmukr(i,2).*0.5.*drho(inter(:,2),2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1),state.h(ic,2).* rhokmukr(i,2).*0.5.*drho(inter(:,1),2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1), state.h(ic,2).*rhokmukr(i,2).*0.5.*drho(inter(:,1),2).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), state.h(ic,2).*rhokmukr(i,2).*0.5.*drho(inter(:,2),2).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ic, state.h(ic,2).*rhokmudkr(i,2).*(dS(ic)).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic, state.h(ic,2).*rhokmudkr(i,2).*(dS(ic)).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic, state.h(ic,2).*rhoskmudkr(i,2).*(dS(ic)).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic, state.h(ic,2).*rhoskmudkr(i,2).*(dS(ic)).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ic,state.dh(ic,2).* rhokmukr(i,2).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic,state.dh(ic,2).* rhokmukr(i,2).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic,state.dh(ic,2).* rhoskmukr(i,2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic,state.dh(ic,2).* rhoskmukr(i,2).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc));


%% m43: dfe/dphi
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,{dfe/dphi},dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}
% derivative of the internal energy due to  the hydration reaction with respect
% to porosity 
d1e=G.cells.volumes./dt.*(state.s(:,1).*state.rho(:,1).*state.u(:,1)+state.s(:,2).*state.rho(:,2).*state.u(:,2)...
    -state.u(:,3).*state.rhoref(:,3));
% the chemical reaction terms (not considered, heatc=0)
d1e=d1e-heatc.*dr_k_phi./0.152.*(state.u(:,3).*0.152-state.u(:,2).*0.044-state.u(:,1).*6.*0.018);
m43=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));



%% m43: dfe/dT
%{
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,{dfe/dT},dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
%}
% derivative of the internal energy due to  the hydration reaction with respect
% to temperature 

d1e=G.cells.volumes.*state.poro./dt.*(state.s(:,1).*drho(:,5).*state.u(:,1)+state.s(:,2).*drho(:,6).*state.u(:,2)...
    +state.s(:,1).*state.rho(:,1).*du(:,3) +state.s(:,2).*state.rho(:,2).*du(:,4))...
    +G.cells.volumes.*(state.poroh)./dt.*(state.rhoref(:,3).*du(:,5))+G.cells.volumes.*(1-rock.poro)./dt.*(state.rhoref(:,5).*du(:,6));
% the chemical reaction terms (not considered, heatc=0)
d1e=d1e-heatc.*dr_k_T./0.152.*(state.u(:,3).*0.152-state.u(:,2).*0.044-state.u(:,1).*6.*0.018);
m44=sparse(1:(ncc),1:(ncc),d1e,(ncc),(ncc));
% advection term
m44=m44...
    +sparse(inter(:,1),ib,state.h(ib,1).* rhokkrdmuT(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* rhokkrdmuT(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.h(ib,1).* rhoskkrdmuT(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* rhoskkrdmuT(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.h(ib,1).* drhokkrmuT(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* drhokkrmuT(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.h(ib,1).* drhoskkrmuT(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.h(ib,1).* drhoskkrmuT(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2),state.h(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,2),5).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1),state.h(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,1),5).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.h(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,1),5).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.h(ib,1).* rhokmukr(i,1).*0.5.*drho(inter(:,2),5).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ib,state.dh(ib,3).* rhokmukr(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.dh(ib,3).* rhokmukr(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),(ncc))...
    +sparse(inter(:,1),ib,state.dh(ib,3).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ib,state.dh(ib,3).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ic,state.h(ic,2).*rhokkrdmuT(i,2).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic, state.h(ic,2).*rhokkrdmuT(i,2).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic, state.h(ic,2).*rhoskkrdmuT(i,2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic,state.h(ic,2).* rhoskkrdmuT(i,2).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ic,state.h(ic,2).* drhokkrmuT(i,2).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic, state.h(ic,2).*drhokkrmuT(i,2).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic,state.h(ic,2).*drhoskkrmuT(i,2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic,state.h(ic,2).* drhoskkrmuT(i,2).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), state.h(ic,2).*rhokmukr(i,2).*0.5.*drho(inter(:,2),6).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1),state.h(ic,2).* rhokmukr(i,2).*0.5.*drho(inter(:,1),6).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),state.h(ic,2).* rhokmukr(i,2).*0.5.*drho(inter(:,1),6).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2),state.h(ic,2).* rhokmukr(i,2).*0.5.*drho(inter(:,2),6).*(-1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),ic,state.dh(ic,4).* rhokmukr(i,2).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),(ncc))...
    +sparse(inter(:,2),ic,state.dh(ic,4).* rhokmukr(i,2).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),(ncc))...
    +sparse(inter(:,1),ic,state.dh(ic,4).* rhoskmukr(i,2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),(ncc))...
    +sparse(inter(:,2),ic,state.dh(ic,4).* rhoskmukr(i,2).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),(ncc))...
    ...
    +sparse(inter(:,1),inter(:,2), KT(i,1).*(-1),(ncc),(ncc))...
    +sparse(inter(:,1),inter(:,1), KT(i,1).*(1),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,2), KT(i,1).*(1),(ncc),(ncc))...
    +sparse(inter(:,2),inter(:,1),KT(i,1).*(-1),(ncc),(ncc));



%% pressure boundary
%{
in the previous Jacobian matrix, we only considered the flux through the
internal face, while the flux through the boundary is not included yet. If
we do not include the flux through the boundary, it means the boundary is
closed.

currently, we only have pressure boundary, on which the gas and liquid
pressures are fixed. The temperature is also fixed on the boundary.
%}
outflux=zeros(1,2);
if (~isempty(bc))
    % pressure on the boundary
    pb=zeros(G.faces.num,2);
    bc.pressure=bc.value;
    pb(bc.face,:)=bc.pressure;
    bb=false(G.faces.num,1);
    bb(bc.face)=true;
    bcell=sum(neighborship(bb,:),2);
    %pb(bc.face,2)=pb(bc.face,1)+(state.pressure(bcell,2)-state.pressure(bcell,1));
    % density and viscosity on the boundary
    [ bc.rho,bc.mu]=rhomu_p_frac_kinetic_h(bc);
    %[drhob,dmub]=D_RHOMU_kinetic_h(bc);
    %specific enthalpy of the boundary
    [bc.h,bc.u,bc.dh,dub]=heat_h(bc,zeros(1,6));
    %density of the boundary
    brho=zeros(G.faces.num,3);
    brho(bc.face,:)=bc.rho;
    %visocity of the boundary
    bmu=zeros(G.faces.num,2);
    bmu(bc.face,:)=bc.mu;
    % specific enthalpy of the boundary
    bh=zeros(G.faces.num,4);
    bh(bc.face,:)=bc.h;
    %relative permeability of the boundary
    bkr=zeros(G.faces.num,2);
    bkr(bc.face,:)=bc.kr;
    % mass fraction of four components in water on the boundary
    bfrac=zeros(G.faces.num,4);
    bfrac(bc.face,:)=bc.frac;
    %find the upwind direction of the boundary face
    upbb=pb(bb,1)-state.pressure(bcell,1)-0.5.*(state.rho(bcell,1)+brho(bb,1)).*(G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'>0;
    upbg=pb(bb,2)-state.pressure(bcell,2)-0.5.*(state.rho(bcell,2)+brho(bb,2)).*(G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'>0;
    upb=zeros(size(upbb));
    upb(upbb)=1;
    upg=zeros(size(upbg));
    upg(upbg)=1;
    % use the relative permeability in the upwind direction
    krb(:,1)=(1-upb).*kr(bcell,1)+upb.*bkr(bb,1);
    krb(:,2)=(1-upg).*kr(bcell,2)+upg.*bkr(bb,2);
    % use the density in the upwind direction
    rhob=brho(bb,:);
    rhob(:,1)=(1-upb).*state.rho(bcell,1)+upb.*brho(bb,1);
    rhob(:,2)=(1-upg).*state.rho(bcell,2)+upg.*brho(bb,2);
    % use the viscosity in the upwind direction
    mub=bmu(bb,:);
    mub(:,1)=(1-upb).*state.mu(bcell,1)+upb.*bmu(bb,1);
    mub(:,2)=(1-upg).*state.mu(bcell,2)+upg.*bmu(bb,2);
    % use the specific enthalpy in the upwind direction
    hb=bh(bb,:);
    hb(:,1)=(1-upb).*state.h(bcell,1)+upb.*bh(bb,1);
    hb(:,2)=(1-upg).*state.h(bcell,2)+upg.*bh(bb,2);
    % use the derivative of the specific enthalpy in the upwind direction
    dhb=state.dh(bcell,:);
    dhb(:,[1 3])=(1-upb).*state.dh(bcell,[1 3]);
    dhb(:,[2 4])=(1-upg).*state.dh(bcell,[2 4]);
    % use the mass fraction of components in water in the upwind direction
    fracb=(1-upb).*state.frac(bcell,:)+upb.*bfrac(bb,:);
    dfracb=(1-upb).*state.dfrac(bcell,:);
    % use the derivative of the density in the upwind direction
    drho(bcell,[1 3 5])=(1-upb).*drho(bcell,[1 3 5]);
    drho(bcell,[2 4 6])=(1-upg).*drho(bcell,[2 4 6]);
    % complete the flux on the boundary face
    state.flux(bc.face,1)=fht(bb).*krb(:,1)./mub(:,1).*(state.pressure(bcell,1)-pb(bb,1)+(0.5.*(state.rho(bcell,1)+brho(bb,1))).*(G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity')./G.faces.areas(bb);
    state.flux(bc.face,2)=fht(bb).*krb(:,2)./mub(:,2).*(state.pressure(bcell,2)-pb(bb,2)+(0.5.*(state.rho(bcell,2)+brho(bb,2))).*(G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity')./G.faces.areas(bb);
    % calculate the total mass flux rate flowing out of the boundary
    outflux(1,1)=(sum(state.flux(bc.face,1).*rhob(:,1).*fracb(:,1).*G.faces.areas(bb))).*dt;
    outflux(1,2)=(sum(state.flux(bc.face,1).*rhob(:,1).*fracb(:,2).*G.faces.areas(bb))+sum(state.flux(bc.face,2).*rhob(:,2).*G.faces.areas(bb))).*dt;
    flux=state.flux;
    state.outflux=outflux;
    % this is similar to previous calculation of the Jacobian matrix due to
    % the flux through the internal face.
    rhokmudkr(bb,1)=(1-upb).*rhob(:,1).*fht(bb).*dkr(bcell,1)./mub(:,1);
    rhokmudkr(bb,2)=(1-upg).*rhob(:,2).*fht(bb).*dkr(bcell,2)./mub(:,2);
    rhoskmudkr(bb,1)=(1-upb).*0.5.*(state.rho(bcell,1)+brho(bb,1)).*rhob(:,1).*fht(bb).*dkr(bcell,1)./mub(:,1);
    rhoskmudkr(bb,2)=(1-upg).*0.5.*(state.rho(bcell,2)+brho(bb,2)).*rhob(:,2).*fht(bb).*dkr(bcell,2)./mub(:,2);
    rhokmukr(bb,1)=rhob(:,1).*fht(bb).*krb(:,1)./mub(:,1);
    rhokmukr(bb,2)=rhob(:,2).*fht(bb).*krb(:,2)./mub(:,2);
    rhokkrdmu(bb,1)=(1-upb).*rhob(:,1).*fht(bb).*krb(:,1)./(state.mu(bcell,1)).^2.*(-dmu(bcell,1));
    rhokkrdmu(bb,2)=(1-upg).*rhob(:,2).*fht(bb).*krb(:,2)./(state.mu(bcell,2)).^2.*(-dmu(bcell,2));
    rhokkrdmuT(bb,1)=(1-upb).*rhob(:,1).*fht(bb).*krb(:,1)./(state.mu(bcell,1)).^2.*(-dmu(bcell,5));
    rhokkrdmuT(bb,2)=(1-upg).*rhob(:,2).*fht(bb).*krb(:,2)./(state.mu(bcell,2)).^2.*(-dmu(bcell,6));
    rhoskkrdmu(bb,1)=(1-upb).*0.5.*(state.rho(bcell,1)+brho(bb,1)).*rhob(:,1).*fht(bb).*krb(:,1)./(state.mu(bcell,1)).^2.*(-dmu(bcell,1));
    rhoskkrdmu(bb,2)=(1-upg).*0.5.*(state.rho(bcell,2)+brho(bb,2)).*rhob(:,2).*fht(bb).*krb(:,2)./(state.mu(bcell,2)).^2.*(-dmu(bcell,2));
    rhoskkrdmuT(bb,1)=(1-upb).*0.5.*(state.rho(bcell,1)+brho(bb,1)).*rhob(:,1).*fht(bb).*krb(:,1)./(state.mu(bcell,1)).^2.*(-dmu(bcell,5));
    rhoskkrdmuT(bb,2)=(1-upg).*0.5.*(state.rho(bcell,2)+brho(bb,2)).*rhob(:,2).*fht(bb).*krb(:,2)./(state.mu(bcell,2)).^2.*(-dmu(bcell,6));

    drhokkrmu(bb,1)=(1-upb).*drho(bcell,1).*fht(bb).*krb(:,1)./(mub(:,1));
    drhokkrmu(bb,2)=(1-upg).*drho(bcell,2).*fht(bb).*krb(:,2)./(mub(:,2));
    drhokkrmuT(bb,1)=(1-upb).*drho(bcell,5).*fht(bb).*krb(:,1)./(mub(:,1));
    drhokkrmuT(bb,2)=(1-upg).*drho(bcell,6).*fht(bb).*krb(:,2)./(mub(:,2));


    drhoskkrmu(bb,1)=0.5.*(1-upb).*(state.rho(bcell,1)+brho(bb,1)).*drho(bcell,1).*fht(bb).*krb(:,1)./(mub(:,1));
    drhoskkrmu(bb,2)=0.5.*(1-upg).*(state.rho(bcell,2)+brho(bb,2)).*drho(bcell,2).*fht(bb).*krb(:,2)./(mub(:,2));
    drhoskkrmuT(bb,1)=0.5.*(1-upb).*(state.rho(bcell,1)+brho(bb,1)).*drho(bcell,5).*fht(bb).*krb(:,1)./(mub(:,1));
    drhoskkrmuT(bb,2)=0.5.*(1-upg).*(state.rho(bcell,2)+brho(bb,2)).*drho(bcell,6).*fht(bb).*krb(:,2)./(mub(:,2));



    rhokkrdmuC(bb,1)=(1-upb).*rhob(:,1).*fht(bb).*krb(:,1)./(mub(:,1)).^2.*(-dmu(bcell,3));
    rhokkrdmuC(bb,2)=(1-upg).*rhob(:,2).*fht(bb).*krb(:,2)./(state.mu(bcell,2)).^2.*(-dmu(bcell,4));
    rhoskkrdmuC(bb,1)=(1-upb).*0.5.*(state.rho(bcell,1)+brho(bb,1)).*rhob(:,1).*fht(bb).*krb(:,1)./(state.mu(bcell,1)).^2.*(-dmu(bcell,3));
    rhoskkrdmuC(bb,2)=(1-upg).*0.5.*(state.rho(bcell,2)+brho(bb,2)).*rhob(:,2).*fht(bb).*krb(:,2)./(state.mu(bcell,2)).^2.*(-dmu(bcell,4));
    drhokkrmuC(bb,1)=(1-upb).*drho(bcell,3).*fht(bb).*krb(:,1)./(mub(:,1));
    drhokkrmuC(bb,2)=(1-upg).*drho(bcell,4).*fht(bb).*krb(:,2)./(mub(:,2));
    drhoskkrmuC(bb,1)=0.5.*(1-upb).*(state.rho(bcell,1)+brho(bb,1)).*drho(bcell,3).*fht(bb).*krb(:,1)./(mub(:,1));
    drhoskkrmuC(bb,2)=0.5.*(1-upg).*(state.rho(bcell,2)+brho(bb,2)).*drho(bcell,4).*fht(bb).*krb(:,2)./(mub(:,2));


    rhoskmukr(bb,1)=0.5.*(state.rho(bcell,1)+brho(bb,1)).*rhob(:,1).*fht(bb).*krb(:,1)./mub(:,1);
    rhoskmukr(bb,2)=0.5.*(state.rho(bcell,2)+brho(bb,2)).*rhob(:,2).*fht(bb).*krb(:,2)./mub(:,2);

%% Jacobian matrix due to bouandry flux
    %this is very similar to the previous calculation for assembling the Jacobian matrix due to the flux through the internal face
   
    m11=m11...
        +sparse(bcell,bcell,fracb(:,1).* rhokkrdmu(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,1).* rhoskkrdmu(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,1).* drhokkrmu(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,1).* drhoskkrmu(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upb).*fracb(:,1).* rhokmudkr(bb,1).*(-dS(bcell)).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upb).*fracb(:,1).* rhoskmudkr(bb,1).*(-dS(bcell)).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,1).* rhokmukr(bb,1).*(1),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,1).* rhokmukr(bb,1).*0.5.*drho(bcell,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc));

    m12=m12...
        +sparse(bcell,bcell,fracb(:,1).* rhokkrdmuC(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...%NEW
        +sparse(bcell,bcell,fracb(:,1).* rhoskkrdmuC(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,1).* drhokkrmuC(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,1).* drhoskkrmuC(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        ...
        +sparse(bcell,bcell,dfracb(:,1).* rhokmukr(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,dfracb(:,1).* rhoskmukr(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upb).*(fracb(:,1)).* rhokmudkr(bb,1).*(dS(bcell)).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upb).*(fracb(:,1)).* rhoskmudkr(bb,1).*(dS(bcell)).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,1).* rhokmukr(bb,1).*0.5.*drho(bcell,3).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc));%NEW

    m14=m14...
        +sparse(bcell,bcell,fracb(:,1).* rhokkrdmuT(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,1).* rhoskkrdmuT(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,1).* drhokkrmuT(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,1).* drhoskkrmuT(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc));

    m21=m21...
        +sparse(bcell,bcell,fracb(:,2).* rhokkrdmu(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,2).* rhoskkrdmu(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,2).* drhokkrmu(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,2).* drhoskkrmu(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upb).*fracb(:,2).* rhokmudkr(bb,1).*(-dS(bcell)).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upb).*fracb(:,2).* rhoskmudkr(bb,1).*(-dS(bcell)).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upg).* rhokmudkr(bb,2).*(-dS(bcell)).*(-pb(bb,2)+state.pressure(bcell,2)),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upg).*rhoskmudkr(bb,2).*(-dS(bcell)).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        ...
        +sparse(bcell,bcell,fracb(:,2).* rhokmukr(bb,1).*(1),(ncc),(ncc))...
        ...%+sparse(bcell,bcell, rhokmukr(bb,2).*(1),(ncc),(ncc))...%pb_g=pb+p_g-p_l;p_g-pb_g=p_l-pb_l
        +sparse(bcell,bcell,fracb(:,2).* rhokmukr(bb,1).*drho(bcell,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc));

    %m22
    m22=m22...
        +sparse(bcell,bcell, rhokkrdmu(bb,2).*(-pb(bb,2)+state.pressure(bcell,2)),(ncc),(ncc))...
        +sparse(bcell,bcell,rhoskkrdmu(bb,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,drhokkrmu(bb,2).*(-pb(bb,2)+state.pressure(bcell,2)),(ncc),(ncc))...
        +sparse(bcell,bcell, drhoskkrmu(bb,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        ...
        +sparse(bcell,bcell,fracb(:,2).* rhokkrdmuC(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...%NEW
        +sparse(bcell,bcell,fracb(:,2).* rhoskkrdmuC(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,2).* drhokkrmuC(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,2).* drhoskkrmuC(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell, (1-upg).*rhokmudkr(bb,2).*(dS(bcell)).*(-pb(bb,2)+state.pressure(bcell,2)),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upg).*rhoskmudkr(bb,2).*(dS(bcell)).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,dfracb(:,2).* rhokmukr(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,dfracb(:,2).* rhoskmukr(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upb).*(fracb(:,2)).* rhokmudkr(bb,1).*(dS(bcell)).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upb).*(fracb(:,2)).* rhoskmudkr(bb,1).*(dS(bcell)).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        ...
        +sparse(bcell,bcell,rhokmukr(bb,2).*(1),(ncc),(ncc))...%pb_g=pb+p_g-p_l;p_g-pb_g=p_l-pb_l
        +sparse(bcell,bcell,rhokmukr(bb,2).*drho(bcell,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,2).* rhokmukr(bb,1).*0.5.*drho(bcell,3).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc));%NEW

    m24=m24...
        +sparse(bcell,bcell,fracb(:,2).* rhokkrdmuT(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,2).* rhoskkrdmuT(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,2).* drhokkrmuT(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,fracb(:,2).* drhoskkrmuT(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        ...
        +sparse(bcell,bcell,rhokkrdmuT(bb,2).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,rhoskkrdmuT(bb,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,drhokkrmuT(bb,2).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell, drhoskkrmuT(bb,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc));

    m41=m41...
        +sparse(bcell,bcell,hb(:,1).* rhokkrdmu(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,1).* rhoskkrdmu(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,1).* drhokkrmu(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,1).* drhoskkrmu(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upb).*hb(:,1).* rhokmudkr(bb,1).*(-dS(bcell)).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upb).*hb(:,1).* rhoskmudkr(bb,1).*(-dS(bcell)).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upg).*hb(:,2).* rhokmudkr(bb,2).*(-dS(bcell)).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upg).*hb(:,2).* rhoskmudkr(bb,2).*(-dS(bcell)).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,1).* rhokmukr(bb,1).*(1),(ncc),(ncc))...
        ...% +sparse(bcell,bcell,hb(:,2).* rhokmukr(bb,2).*(1),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,1).* rhokmukr(bb,1).*0.5.*drho(bcell,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc));

    m42=m42...
        +sparse(bcell,bcell,hb(:,2).* rhokkrdmu(bb,2).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,2).* rhoskkrdmu(bb,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,2).* drhokkrmu(bb,2).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,2).* drhoskkrmu(bb,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upg).*hb(:,2).* rhokmudkr(bb,2).*(dS(bcell)).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upg).*hb(:,2).* rhoskmudkr(bb,2).*(dS(bcell)).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upb).*hb(:,1).* rhokmudkr(bb,1).*(dS(bcell)).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,(1-upb).*hb(:,1).* rhoskmudkr(bb,1).*(dS(bcell)).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        ...
        +sparse(bcell,bcell,hb(:,2).* rhokmukr(bb,2).*drho(bcell,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,dhb(:,2).* rhokmukr(bb,2).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,dhb(:,2).* rhoskmukr(bb,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,2).* rhokmukr(bb,2).*(1),(ncc),(ncc));

    m44=m44...
        +sparse(bcell,bcell,hb(:,1).* rhokkrdmuT(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,1).* rhoskkrdmuT(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,1).* drhokkrmuT(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,1).* drhoskkrmuT(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        ...
        +sparse(bcell,bcell,hb(:,2).* rhokkrdmuT(bb,2).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,2).* rhoskkrdmuT(bb,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,2).* drhokkrmuT(bb,2).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,hb(:,2).* drhoskkrmuT(bb,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        ...
        +sparse(bcell,bcell,dhb(:,3).* rhokmukr(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,dhb(:,3).* rhoskmukr(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc))...
        ...
        +sparse(bcell,bcell,dhb(:,4).* rhokmukr(bb,2).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),(ncc))...
        +sparse(bcell,bcell,dhb(:,4).* rhoskmukr(bb,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),(ncc));



end



%due to wells
nw=length(W);
wsat=state.s;

for kw=1:nw
    wc=W(kw).cells;
    wsat(wc,:)=repmat(W(kw).compi,[length(wc),1]);
end

krw=fluid.relperm(wsat);%%%%%%be careful about the the relative permeability;state.s


wCol1=sparse([],[],[],ncc,nw);
wCol2=sparse([],[],[],ncc,nw);
wCol3=sparse([],[],[],ncc,nw);
wCol4=sparse([],[],[],ncc,nw);
wRow1=sparse([],[],[],nw,ncc);
wRow2=sparse([],[],[],nw,ncc);
wRow3=sparse([],[],[],nw,ncc);
wRow4=sparse([],[],[],nw,ncc);
wCor=sparse([],[],[],nw,nw);


h_W=zeros(nw,4);
dh_W=zeros(nw,6);
for ki= 1 : nw



    if W(ki).compi(1)==1 %inject brine

        W(ki).dfrac=[0 0 0 0];

        W(ki).rhoref=state.rhoref;

        [W(ki).rho,W(ki).mu]=rhomu_p_frac_kinetic_h(W(ki));

        [drhow,dmuw]=D_RHOMU_kinetic_h(W(ki));
        [h_W(ki,:),~,dh_W(ki,:),~]=heat_h(W(ki),drhow);
        wc=W(ki).cells;




        %for water component
        WG=-( W(ki).frac(1).*drhow(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)...
            +W(ki).frac(1).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1).^2.*(-dmuw(1)))...
            .*(W(ki).pressure(:,1)-state.pressure(wc,1)+W(ki).rho(1).*W(ki).dZ.*norm(gravity))...
            -W(ki).frac(1).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)...
            -W(ki).frac(1).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1).*drhow(1).*W(ki).dZ.*norm(gravity);
        FG=W(ki).frac(1).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1);
        Wdg=sum(WG);
        wCol1=wCol1+sparse(wc,ki,WG,(ncc),(nw));
        wRow1=wRow1+sparse(ki,wc,-FG,(nw),(ncc));
        wCor=wCor+sparse(ki,ki,-Wdg,(nw),(nw));

        m11=m11+sparse(wc,wc,FG,(ncc),(ncc));

        %for co2 component
        WG=-(W(ki).frac(2).*drhow(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)...
            +W(ki).frac(2).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1).^2.*(-dmuw(1)))...
            .*(W(ki).pressure(:,1)-state.pressure(wc,1)+W(ki).rho(1).*W(ki).dZ.*norm(gravity))...
            -W(ki).frac(2).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)...
            -W(ki).frac(2).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1).*drhow(1).*W(ki).dZ.*norm(gravity);
        FG=W(ki).frac(2).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1);
        Wdg=sum(WG);
        wCol2=wCol2+sparse(wc,ki,WG,(ncc),(nw));
        wRow1=wRow1+sparse(ki,wc,-FG,(nw),(ncc));
        wCor=wCor+sparse(ki,ki,-Wdg,(nw),(nw));

        m21=m21+sparse(wc,wc,FG,(ncc),(ncc));

        %for mineral component
        WG=-(W(ki).frac(3).*drhow(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)...
            +W(ki).frac(3).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1).^2.*(-dmuw(1)))...
            .*(W(ki).pressure(:,1)-state.pressure(wc,1)+W(ki).rho(1).*W(ki).dZ.*norm(gravity))...
            -W(ki).frac(3).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)...
            -W(ki).frac(3).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1).*drhow(1).*W(ki).dZ.*norm(gravity);
        FG=W(ki).frac(3).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1);
        Wdg=sum(WG);

        wRow1=wRow1+sparse(ki,wc,-FG,(nw),(ncc));
        wCor=wCor+sparse(ki,ki,-Wdg,(nw),(nw));


        %for heat component
        WG=-(dh_W(ki,1).*rhow(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)+h_W(ki,1).*drhow(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)...
            +h_W(ki,1).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1).^2.*(-dmuw(1)))...
            .*(W(ki,1).pressure(:,1)-state.pressure(wc,1)+W(ki).rho(1).*W(ki).dZ.*norm(gravity))...
            -h_W(ki,1).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)...
            -h_W(ki,1).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1).*drhow(1).*W(ki).dZ.*norm(gravity);
        FG=h_W(ki,1).*W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1);

        wCol4=wCol4+sparse(wc,ki,WG,(ncc),(nw));



        m41=m41+sparse(wc,wc,FG,(ncc),(ncc));

    else

       wc=W(ki).cells;
        W(ki).frac=[0 0 0 1];
        W(ki).dfrac=[0 0 0 0];
        W(ki).rhoref=state.rhoref;
        Wrep=W(ki);
        Wrep.pressure=[max(state.pressure0(wc,1)),max(state.pressure0(wc,1))];
        [W(ki).rho,W(ki).mu]=rhomu_p_frac_kinetic_h(Wrep);
        Wrep.rho=W(ki).rho;
         Wrep.mu=W(ki).mu;
        Wrep.dsp=0;
        Wrep.dcp=0;
        [drhow,dmuw]=D_RHOMU_kinetic_h(Wrep);
         drhow=zeros(size(drhow));
         dmuw=zeros(size(dmuw));
         
 
        [h_W(ki,:),~,dh_W(ki,:),~]=heat_h(W(ki),drhow);
       
        %for co2 species
        WG=-(W(ki).frac(4).*drhow(2).*W(ki).WI.*krw(wc,2)./W(ki).mu(2)...
            +W(ki).frac(4).*W(ki).rho(2).*W(ki).WI.*krw(wc,2)./W(ki).mu(2).^2.*(-dmuw(2)))...
            .*(W(ki).pressure(:,2)-state.pressure(wc,2)+W(ki).rho(2).*W(ki).dZ.*norm(gravity))...
            -W(ki).frac(4).*W(ki).rho(2).*W(ki).WI.*krw(wc,2)./W(ki).mu(2)...
            -W(ki).frac(4).*W(ki).rho(2).*W(ki).WI.*krw(wc,2)./W(ki).mu(2).*drhow(2).*W(ki).dZ.*norm(gravity);
        FG=W(ki).frac(4).*W(ki).rho(2).*W(ki).WI.*krw(wc,2)./W(ki).mu(2);
        Wdg=sum(WG);
        wCol2=wCol2+sparse(wc,ki,WG,(ncc),(nw));
        wRow2=wRow2+sparse(ki,wc,-FG,(nw),(ncc));
        wCor=wCor+sparse(ki,ki,-Wdg,(nw),(nw));
        m22=m22+sparse(wc,wc,FG,(ncc),(ncc));

        %for heat component
        WG=-(dh_W(ki,2).*W(ki).rho(2).*W(ki).WI.*krw(wc,2)./W(ki).mu(2)+h_W(ki,2).*drhow(2).*W(ki).WI.*krw(wc,2)./W(ki).mu(2)...
            +h_W(ki,2).*W(ki).rho(2).*W(ki).WI.*krw(wc,2)./W(ki).mu(2).^2.*(-dmuw(2)))...
            .*(W(ki).pressure(:,2)-state.pressure(wc,2)+W(ki).rho(2).*W(ki).dZ.*norm(gravity))...
            -h_W(ki,2).*W(ki).rho(2).*W(ki).WI.*krw(wc,2)./W(ki).mu(2)...
            -h_W(ki,2).*W(ki).rho(2).*W(ki).WI.*krw(wc,2)./W(ki).mu(2).*drhow(2).*W(ki).dZ.*norm(gravity);
        FG=h_W(ki,2).*W(ki).rho(2).*W(ki).WI.*krw(wc,2)./W(ki).mu(2);
        wCol4=wCol4+sparse(wc,ki,WG,(ncc),(nw));
        m42=m42+sparse(wc,wc,FG,(ncc),(ncc));


    end
end








%% RHS: right hand side of the Neton equation L*x=R; i.e., residual vector
% R=[-fw;-fc;-fe;-fW]=[R1;R2;R3;R4]
accum=zeros(length(G.cells.volumes),4);
re1=1;
re2=1-re1;
d2e=-G.cells.volumes./dt.*...
    (state.rho(:,1).*state.s(:,1).*state.frac(:,1).*state.poro-state.rho0(:,1).*state.s0(:,1).*state.frac0(:,1).*state.poro0...
    +state.rhoref(:,3).*(state.poroh).*state.omega_t(1,6)-state.rhoref(:,3).*(state.poroh0).*state.omega_t(1,6));

d1e=-(G.cells.volumes./dt.*...
    (state.rho(:,1).*state.s(:,1).*state.frac(:,1).*state.poro-state.rho0(:,1).*state.s0(:,1).*state.frac0(:,1).*state.poro0...
    )+1.*6.*18./152.*r_k);
d1e=(re1.*d1e+re2.*d2e);

accum(:,1)=-G.cells.volumes./dt.*...
    (state.rho(:,1).*state.s(:,1).*state.frac(:,1).*state.poro...
    +state.rhoref(:,3).*(state.poroh).*state.omega_t(1,6));
R1=sparse(1:(ncc),1,d1e,(ncc),1);
R1=R1...
    +sparse(inter(:,1),1,-state.frac(ib,1).* rhokmukr(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),1)...
    +sparse(inter(:,2),1,-state.frac(ib,1).* rhokmukr(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),1)...
    +sparse(inter(:,1),1,-state.frac(ib,1).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),1)...
    +sparse(inter(:,2),1,-state.frac(ib,1).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),1)...
    ...
    +sparse(inter(:,1),1, sphirhod(i,1).*(state.frac(inter(:,2),1)-state.frac(inter(:,1),1)),(ncc),1)...
    +sparse(inter(:,2),1, sphirhod(i,1).*(state.frac(inter(:,1),1)-state.frac(inter(:,2),1)),(ncc),1);


d2e=-G.cells.volumes./dt.*...
    (state.rho(:,1).*state.s(:,1).*state.frac(:,2).*state.poro-state.rho0(:,1).*state.s0(:,1).*state.frac0(:,2).*state.poro0...
    +state.rho(:,2).*state.s(:,2).*state.poro-state.rho0(:,2).*state.s0(:,2).*state.poro0...
    +state.rhoref(:,3).*(state.poroh).*state.omega_t(2,6)-state.rhoref(:,3).*(state.poroh0).*state.omega_t(2,6));

d1e=-(G.cells.volumes./dt.*...
    (state.rho(:,1).*state.s(:,1).*state.frac(:,2).*state.poro-state.rho0(:,1).*state.s0(:,1).*state.frac0(:,2).*state.poro0...
    +state.rho(:,2).*state.s(:,2).*state.poro-state.rho0(:,2).*state.s0(:,2).*state.poro0...
    ) +1.*44./152.*r_k);
d1e=(re1.*d1e+re2.*d2e);

accum(:,2)=-G.cells.volumes./dt.*...
    (state.rho(:,1).*state.s(:,1).*state.frac(:,2).*state.poro...
    +state.rho(:,2).*state.s(:,2).*state.poro...
    +state.rhoref(:,3).*(state.poroh).*state.omega_t(2,6));
R2=sparse(1:(ncc),1,d1e,(ncc),1);

R2=R2...
    +sparse(inter(:,1),1,-(state.frac(ib,2)).* rhokmukr(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),1)...
    +sparse(inter(:,2),1,-(state.frac(ib,2)).* rhokmukr(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),1)...
    +sparse(inter(:,1),1,-(state.frac(ib,2)).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),1)...
    +sparse(inter(:,2),1,-(state.frac(ib,2)).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),1)...
    ...
    +sparse(inter(:,1),1,-rhokmukr(i,2).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),1)...
    +sparse(inter(:,2),1,- rhokmukr(i,2).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),1)...
    +sparse(inter(:,1),1,-rhoskmukr(i,2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),1)...
    +sparse(inter(:,2),1,-rhoskmukr(i,2).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),1)...
    ...
    +sparse(inter(:,1),1, sphirhod(i,1).*(state.frac(inter(:,2),2)-state.frac(inter(:,1),2)),(ncc),1)...
    +sparse(inter(:,2),1, sphirhod(i,1).*(state.frac(inter(:,1),2)-state.frac(inter(:,2),2)),(ncc),1);




d2e=-(G.cells.volumes./dt.*...
    (state.rhoref(3).*(1-state.poro-(1-rock.poro))-state.rhoref(3).*(1-state.poro0-(1-rock.poro)))-...
    (state.species(:,6)-state.species0(:,6))./dt);
d1e=-(G.cells.volumes./dt.*...
    (state.rhoref(3).*(1-state.poro-(1-rock.poro))-state.rhoref(3).*(1-state.poro0-(1-rock.poro)))-r_k);
d1e=(re1.*d1e+re2.*d2e);

accum(:,3)=-(G.cells.volumes./dt.*...
    (state.rhoref(3).*(1-state.poro-(1-rock.poro)))-...
    (state.species(:,6))./dt);
R3=sparse(1:(ncc),1,d1e,(ncc),1);


d1e=-G.cells.volumes./dt.*...
    (state.rho(:,1).*state.s(:,1).*state.u(:,1).*state.poro-state.rho0(:,1).*state.s0(:,1).*state.u0(:,1).*state.poro0...
    +state.rho(:,2).*state.s(:,2).*state.poro.*state.u(:,2)-state.rho0(:,2).*state.s0(:,2).*state.poro0.*state.u0(:,2)...
    +state.rhoref(:,3).*(state.poroh).*state.u(:,3)-state.rhoref(:,3).*(state.poroh0).*state.u0(:,3)...
    +state.rhoref(:,5).*(1-rock.poro).*state.u(:,4)-state.rhoref(:,5).*(1-rock.poro).*state.u0(:,4));
d1e=d1e+heatc.*r_k./0.152.*(state.u(:,3).*0.152-state.u(:,2).*0.044-state.u(:,1).*6.*0.018);


accum(:,4)=-G.cells.volumes./dt.*...
    (state.rho(:,1).*state.s(:,1).*state.u(:,1).*state.poro...
    +state.rho(:,2).*state.s(:,2).*state.poro.*state.u(:,2)...
    +state.rhoref(:,3).*(state.poroh).*state.u(:,3)...
    +state.rhoref(:,5).*(1-rock.poro).*state.u(:,4));

R4=sparse(1:(ncc),1,d1e,(ncc),1);

R4=R4...
    +sparse(inter(:,1),1,-(state.h(ib,1)).* rhokmukr(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),1)...
    +sparse(inter(:,2),1,-(state.h(ib,1)).* rhokmukr(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),1)...
    +sparse(inter(:,1),1,-(state.h(ib,1)).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),1)...
    +sparse(inter(:,2),1,-(state.h(ib,1)).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),1)...
    ...
    +sparse(inter(:,1),1,-(state.h(ic,2)).*rhokmukr(i,2).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),1)...
    +sparse(inter(:,2),1,- (state.h(ic,2)).*rhokmukr(i,2).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),1)...
    +sparse(inter(:,1),1,-(state.h(ic,2)).*rhoskmukr(i,2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),1)...
    +sparse(inter(:,2),1,-(state.h(ic,2)).*rhoskmukr(i,2).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),1)...
    ...
    +sparse(inter(:,1),1, KT(i,1).*(state.Tk(inter(:,2),1)-state.Tk(inter(:,1),1)),(ncc),1)...
    +sparse(inter(:,2),1, KT(i,1).*(state.Tk(inter(:,1),1)-state.Tk(inter(:,2),1)),(ncc),1);




%% boundary effect on RHS (right hand side of the Neton equation L*x=R; i.e., residual vector)


if ~isempty(bc)&&sum( bc.tb==1)

    R1=R1...
        +sparse(bcell,1,-fracb(:,1).* rhokmukr(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),1)...
        +sparse(bcell,1,-fracb(:,1).* rhoskmukr(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),1);



    R2=R2...
        +sparse(bcell,1,-(fracb(:,2)).* rhokmukr(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),1)...
        +sparse(bcell,1,-(fracb(:,2)).* rhoskmukr(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),1)...
        +sparse(bcell,1,-rhokmukr(bb,2).*(-pb(bb,2)+state.pressure(bcell,2)),(ncc),1)...
        +sparse(bcell,1,-rhoskmukr(bb,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),1);


    R4=R4...
        +sparse(bcell,1,-(hb(:,1)).* rhokmukr(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),1)...
        +sparse(bcell,1,-(hb(:,1)).* rhoskmukr(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),1)...
        +sparse(bcell,1,-hb(:,2).*rhokmukr(bb,2).*(-pb(bb,2)+state.pressure(bcell,2)),(ncc),1)...
        +sparse(bcell,1,-hb(:,2).*rhoskmukr(bb,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),1);

end


%% effect of injection well on RHS (right hand side of the Neton equation L*x=R; i.e., residual vector)

rW=sparse([],[],[],(nw),1);
for ki= 1 : nw

    wc=W(ki).cells;
    if W(ki).compi(1)==1 %inject brine


        % for co2 comonent in the grid
        RF=(W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)).*(W(ki).pressure(:,1)-state.pressure(wc,1)+W(ki).rho(1).*W(ki).dZ.*norm(gravity));
        R2=R2+sparse(wc,1,W(ki).frac(2).*RF,(ncc),1);
        % for water component in the grid
        R1=R1+sparse(wc,1,W(ki).frac(1).*RF,(ncc),1);
        % for mineral component in the grid

        R4=R4+sparse(wc,1,h_W(ki,1).*RF,(ncc),1);
        %for injection well  both species
        RW=W(ki).val.*state.rhoref(1)-sum(RF);
        rW=rW+sparse(ki,1,RW,(nw),1);



    else
        % for co2 species in the grid
        RF=(W(ki).rho(2).*W(ki).WI.*krw(wc,2)./W(ki).mu(2)).*(W(ki).pressure(:,2)-state.pressure(wc,2)+W(ki).rho(2).*W(ki).dZ.*norm(gravity));
        R2=R2+sparse(wc,1,W(ki).frac(4).*RF,(ncc),1);
        R4=R4+sparse(wc,1,h_W(ki,2).*RF,(ncc),1);
        %for injection well  both species
        RW=W(ki).val.*state.rhoref(2)-sum(RF);
        rW=rW+sparse(ki,1,RW,(nw),1);
    end
end


%{ 
 jacobian matrix L=[dfw/dpl,dfw/dpg,dfw/dphi,dfw/dT,dfw/dpbh;
                  dfc/dpl,dfc/dpg,dfc/dphi,dfc/dT,dfc/dpbh;
                  dfHYD/dpl,dfHYD/dpg,dfHYD/dphi,dfHYD/dT,dfHYD/dpbh;
                  dfe/dpl,dfe/dpg,dfe/dphi,dfe/dT,dfe/dpbh;  
                  dfW/dpl,dfW/dpg,dfW/dphi,dfW/dT,dfW/dpbh]
 residual vector R=[-fw;-fc,-fe,-fW]
here, fw=governing equation of water component
      fc=governing equation of CO2 component
      fHYD=governing equation of hydrate component
      fe=governing equation of energy
      fW=governing equation of CO2 mass in the injection well
      pl=liquid pressure
      pg=gas pressure
      phi=volume fraction of the fluid phase
      T=temperature
      pbh=botttom hole pressure
%}
L=[m11,m12,m13,m14,wCol1;...
    m21,m22,m23,m24,wCol2;...
    m31,m32,m33,m34,wCol3;...
    m41,m42,m43,m44,wCol4;...
    wRow1,wRow2,wRow3,wRow4,wCor];
R=[R1;R2;R3;R4;rW];

%% update the mass of the charge component, Cl component, CO2 component and water component
COMPONENT=zeros(size(state.component));
%% update the mass of the charge component
% the initial mass of the charge component (this mass is divided by dt, later it will be multiplied by dt)
deo=G.cells.volumes./dt.*...
    (state.rho0(:,1).*state.s0(:,1).*state.frac0(:,4).*state.poro0);
Ro=sparse(1:(ncc),1,deo,(ncc),1);
% mass change due flux through  interfaces
Ro=Ro+sparse(inter(:,1),1,-state.frac(ib,4).* rhokmukr(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),1)...
    +sparse(inter(:,2),1,-state.frac(ib,4).* rhokmukr(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),1)...
    +sparse(inter(:,1),1,-state.frac(ib,4).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),1)...
    +sparse(inter(:,2),1,-state.frac(ib,4).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),1)...
    ...
    +sparse(inter(:,1),1, sphirhod(i,1).*(state.frac(inter(:,2),4)-state.frac(inter(:,1),4)),(ncc),1)...
    +sparse(inter(:,2),1, sphirhod(i,1).*(state.frac(inter(:,1),4)-state.frac(inter(:,2),4)),(ncc),1);

% mass change due flux through boundary  interfaces
if ~isempty(bc)

    Ro=Ro...
        +sparse(bcell,1,-fracb(:,4).* rhokmukr(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),1)...
        +sparse(bcell,1,-fracb(:,4).* rhoskmukr(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),1);
end
% mass change due to the injection well
for ki= 1 : nw

    wc=W(ki).cells;
    if W(ki).compi(1)==1 %inject brine

        RF=(W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)).*(W(ki).pressure(:,1)-state.pressure(wc,1)+W(ki).rho(1).*W(ki).dZ.*norm(gravity));
        Ro=Ro+sparse(wc,1,W(ki).frac(4).*RF,(ncc),1);
    end
end

COMPONENT(:,4)=Ro.*dt;
%% update the mass of the Cl component
% initial mass
deo=G.cells.volumes./dt.*...
    (state.rho0(:,1).*state.s0(:,1).*state.frac0(:,3).*state.poro0+state.rhoref(:,4).*(state.porona0));
Ro=sparse(1:(ncc),1,deo,(ncc),1);
% mass change due to flux through interfaces
Ro=Ro+sparse(inter(:,1),1,-state.frac(ib,3).* rhokmukr(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),1)...
    +sparse(inter(:,2),1,-state.frac(ib,3).* rhokmukr(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),1)...
    +sparse(inter(:,1),1,-state.frac(ib,3).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),1)...
    +sparse(inter(:,2),1,-state.frac(ib,3).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),1)...
    ...
    +sparse(inter(:,1),1, sphirhod(i,1).*(state.frac(inter(:,2),3)-state.frac(inter(:,1),3)),(ncc),1)...
    +sparse(inter(:,2),1, sphirhod(i,1).*(state.frac(inter(:,1),3)-state.frac(inter(:,2),3)),(ncc),1);
% mass change due to boundary faces
if ~isempty(bc)

    Ro=Ro...
        +sparse(bcell,1,-fracb(:,3).* rhokmukr(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),1)...
        +sparse(bcell,1,-fracb(:,3).* rhoskmukr(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),1);
end
% mass change due to injection well
for ki= 1 : nw
    wc=W(ki).cells;
    if W(ki).compi(1)==1 %inject brine
        RF=(W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)).*(W(ki).pressure(:,1)-state.pressure(wc,1)+W(ki).rho(1).*W(ki).dZ.*norm(gravity));
        Ro=Ro+sparse(wc,1,W(ki).frac(3).*RF,(ncc),1);
    end
end

COMPONENT(:,3)=Ro.*dt;
%% update the mass of the CO2 component
% initial mass
deo=G.cells.volumes./dt.*...
    (state.rho0(:,1).*state.s0(:,1).*state.frac0(:,2).*state.poro0+state.rho0(:,2).*state.s0(:,2).*state.poro0...
    +state.rhoref(:,3).*state.poroh0.*state.omega_t(2,6));
Ro=sparse(1:(ncc),1,deo,(ncc),1);
% mass change due to flux through interfaces
Ro=Ro+sparse(inter(:,1),1,-(state.frac(ib,2)).* rhokmukr(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),1)...
    +sparse(inter(:,2),1,-(state.frac(ib,2)).* rhokmukr(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),1)...
    +sparse(inter(:,1),1,-(state.frac(ib,2)).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),1)...
    +sparse(inter(:,2),1,-(state.frac(ib,2)).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),1)...
    ...
    +sparse(inter(:,1),1,-rhokmukr(i,2).*(state.pressure(inter(:,1),2)-state.pressure(inter(:,2),2)),(ncc),1)...
    +sparse(inter(:,2),1,- rhokmukr(i,2).*(state.pressure(inter(:,2),2)-state.pressure(inter(:,1),2)),(ncc),1)...
    +sparse(inter(:,1),1,-rhoskmukr(i,2).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),1)...
    +sparse(inter(:,2),1,-rhoskmukr(i,2).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),1)...
    ...
    +sparse(inter(:,1),1, sphirhod(i,1).*(state.frac(inter(:,2),2)-state.frac(inter(:,1),2)),(ncc),1)...
    +sparse(inter(:,2),1, sphirhod(i,1).*(state.frac(inter(:,1),2)-state.frac(inter(:,2),2)),(ncc),1);

% mass change due to boundary faces
if ~isempty(bc)

    Ro=Ro...
        +sparse(bcell,1,-fracb(:,2).* rhokmukr(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),1)...
        +sparse(bcell,1,-fracb(:,2).* rhoskmukr(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),1)...
        +sparse(bcell,1,-rhokmukr(bb,2).*(-pb(bb,2)+state.pressure(bcell,2)),(ncc),1)...
        +sparse(bcell,1,-rhoskmukr(bb,2).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),1);
end
% mass change due to the injection well
for ki= 1 : nw

    wc=W(ki).cells;
    if W(ki).compi(1)==1 %inject brine

        RF=(W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)).*(W(ki).pressure(:,1)-state.pressure(wc,1)+W(ki).rho(1).*W(ki).dZ.*norm(gravity));
        Ro=Ro+sparse(wc,1,W(ki).frac(2).*RF,(ncc),1);

    else
        % for co2 species in the grid
        RF=(W(ki).rho(2).*W(ki).WI.*krw(wc,2)./W(ki).mu(2)).*(W(ki).pressure(:,2)-state.pressure(wc,2)+W(ki).rho(2).*W(ki).dZ.*norm(gravity));
        Ro=Ro+sparse(wc,1,W(ki).frac(4).*RF,(ncc),1);

    end
end

COMPONENT(:,2)=Ro.*dt;

%% update the mass of the water component
% initial mass
deo=G.cells.volumes./dt.*...
    (state.rho0(:,1).*state.s0(:,1).*state.frac0(:,1).*state.poro0+state.rhoref(:,3).*state.poroh0.*state.omega_t(1,6));
Ro=sparse(1:(ncc),1,deo,(ncc),1);
% mass change due to flux through interfaces
Ro=Ro+sparse(inter(:,1),1,-state.frac(ib,1).* rhokmukr(i,1).*(state.pressure(inter(:,1),1)-state.pressure(inter(:,2),1)),(ncc),1)...
    +sparse(inter(:,2),1,-state.frac(ib,1).* rhokmukr(i,1).*(state.pressure(inter(:,2),1)-state.pressure(inter(:,1),1)),(ncc),1)...
    +sparse(inter(:,1),1,-state.frac(ib,1).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,2),:)-G.cells.centroids(inter(:,1),:))*gravity'),(ncc),1)...
    +sparse(inter(:,2),1,-state.frac(ib,1).* rhoskmukr(i,1).*((G.cells.centroids(inter(:,1),:)-G.cells.centroids(inter(:,2),:))*gravity'),(ncc),1)...
    ...
    +sparse(inter(:,1),1, sphirhod(i,1).*(state.frac(inter(:,2),1)-state.frac(inter(:,1),1)),(ncc),1)...
    +sparse(inter(:,2),1, sphirhod(i,1).*(state.frac(inter(:,1),1)-state.frac(inter(:,2),1)),(ncc),1);


% mass change due to boundary faces
if ~isempty(bc)
    Ro=Ro...
        +sparse(bcell,1,-fracb(:,1).* rhokmukr(bb,1).*(-pb(bb,1)+state.pressure(bcell,1)),(ncc),1)...
        +sparse(bcell,1,-fracb(:,1).* rhoskmukr(bb,1).*((G.faces.centroids(bb,:)-G.cells.centroids(bcell,:))*gravity'),(ncc),1);
end
% mass change due to the injection well
for ki= 1 : nw

    wc=W(ki).cells;
    if W(ki).compi(1)==1 %inject brine

        RF=(W(ki).rho(1).*W(ki).WI.*krw(wc,1)./W(ki).mu(1)).*(W(ki).pressure(:,1)-state.pressure(wc,1)+W(ki).rho(1).*W(ki).dZ.*norm(gravity));
        Ro=Ro+sparse(wc,1,W(ki).frac(1).*RF,(ncc),1);
    end
end
COMPONENT(:,1)=Ro.*dt;
end




function [ib,ic]=upwind(G,state)
% get the upwind cell of each face
neighbors = getNeighbourship(G, 'Topological', true);
intern = all(neighbors ~= 0, 2);

gb=(state.pressure(neighbors(intern,2),1)-state.pressure(neighbors(intern,1),1))...
    -(state.rho(neighbors(intern,1),1)+state.rho(neighbors(intern,2),1))./2.*((G.cells.centroids(neighbors( intern,2),:)-G.cells.centroids(neighbors( intern,1),:))*gravity');
gc=(state.pressure(neighbors(intern,2),2)-state.pressure(neighbors(intern,1),2))...
    -(state.rho(neighbors(intern,1),2)+state.rho(neighbors(intern,2),2))./2.*((G.cells.centroids(neighbors( intern,2),:)-G.cells.centroids(neighbors( intern,1),:))*gravity');
% check the water flux direction
c = gb>0;
N       = neighbors(intern,:);
N(c, :) = N(c, [2,1]);
% upwind cell  for water flux
ib= N(:, 1);
% check the CO2 flux direction
c = gc>0;
N       = neighbors(intern,:);
N(c, :) = N(c, [2,1]);
% Upwind cell for CO2 flux.
ic= N(:, 1);

end

function [kb,kc,KT]=DISP_H(state,G,rock)
%{
The detailed description can be found in pages 131 -134 of LIE, Knut-Andreas. An introduction to reservoir simulation using MATLAB/GNU Octave: User guide for the MATLAB Reservoir Simulation Toolbox (MRST). Cambridge University Press, 2019.
--------------
here we calculate 
kb: the half mass transmissivity due to dispersion  in brine of each cell; if this cell has 6 faces, then kb contains 6 values for 6 faces, respectively. 
kc: the half mass transmissivity due to dispersion in gas of each cell; if this cell has 6 faces, then kc contains 6 values for 6 faces, respectively. 
kT: the half heat transmissivity of each cell; if this cell has 6 faces, then kT contains 6 values for 6 faces, respectively. 

We only give a demonstrating example of how to derive kb of cell i with respect face f (interface between cells i and j). 

kb_i=A_ij*D_i*c_if/|c_if|^2*n_f, 
where: subscripts i and j denote two ajacent cells,
         A_ij= face area of interface between cell i and j;D_i =[D_i_xx,D_i_xy,D_i_xz;D_i_yx,D_i_yy,D_i_yz;D_i_zx,D_i_zy,D_i_zz] is the dispersion tensor
          c_if= the vector from the center of cell i to the center of the interface between cells i and j, 
       and n_f =the unit vector representing the interface direction pointing out of cell i. 

The flux due to dispersion is calculated as 

F_ij=s_l*phi*rho_l*1/(1/kb_i+1/kb_j)*(X_i-X_j)
here,   s_l= water saturation; phi=porosity; rho_l=water density; X= concentration 
%}
[neighborship, ~] = getNeighbourship(G, 'Topological', true);

i  = all(neighborship ~= 0, 2);
% water discharge (q_l) of each face
flux_face_vector_brine=state.flux(:,1)./G.faces.areas(:).*G.faces.normals;
% CO2 discharge (q_g) of each face
flux_face_vector_gas=state.flux(:,2)./G.faces.areas(:).*G.faces.normals;
% water velocity(v_l) of each cell in x, y and z directions(v_l=q_l/phi)
flux_cell_brine(:,1)=1./rock.poro.*accumarray([neighborship(i,1); neighborship(i,2)],[ flux_face_vector_brine(i,1);flux_face_vector_brine(i,1)])./accumarray([neighborship(i,1); neighborship(i,2)], 1);
flux_cell_brine(:,2)=1./rock.poro.*accumarray([neighborship(i,1); neighborship(i,2)],[ flux_face_vector_brine(i,2);flux_face_vector_brine(i,2)])./accumarray([neighborship(i,1); neighborship(i,2)], 1);
flux_cell_brine(:,3)=1./rock.poro.*accumarray([neighborship(i,1); neighborship(i,2)],[ flux_face_vector_brine(i,3);flux_face_vector_brine(i,3)])./accumarray([neighborship(i,1); neighborship(i,2)], 1);
% CO2 velocity (v_g) of each cell in x, y and z directions(v_g=q_g/phi)
flux_cell_gas(:,1)=1./rock.poro.*accumarray([neighborship(i,1); neighborship(i,2)],[ flux_face_vector_gas(i,1);flux_face_vector_gas(i,1)])./accumarray([neighborship(i,1); neighborship(i,2)], 1);
flux_cell_gas(:,2)=1./rock.poro.*accumarray([neighborship(i,1); neighborship(i,2)],[ flux_face_vector_gas(i,2);flux_face_vector_gas(i,2)])./accumarray([neighborship(i,1); neighborship(i,2)], 1);
flux_cell_gas(:,3)=1./rock.poro.*accumarray([neighborship(i,1); neighborship(i,2)],[ flux_face_vector_gas(i,3);flux_face_vector_gas(i,3)])./accumarray([neighborship(i,1); neighborship(i,2)], 1);
% absolute water velocity
vabs(:,1)=sqrt(sum(flux_cell_brine.^2,2));
% absolute CO2 velocity
vabs(:,2)=sqrt(sum(flux_cell_gas.^2,2));
%% half mass transmissivity due to dispersion in water phase
al=state.al(1); %longitudinal dispersivity
at=state.at(1); %lateral dispersivity
dm=state.dm(1);
taubc=rock.poro.^(4/3).*(state.s0(:,1)).^(10/3);
tau=taubc;
% dispersion tensor of the water phase
Kb(:,1)=al.*vabs(:,1)+(al-at).*( flux_cell_brine(:,1)).^2./vabs(:,1)+tau.*dm;%xx
Kb(:,2)=(al-at).* flux_cell_brine(:,1).* flux_cell_brine(:,2)./vabs(:,1);%xy
Kb(:,3)=(al-at).* flux_cell_brine(:,1).* flux_cell_brine(:,3)./vabs(:,1);%xz
Kb(:,4)=al.*vabs(:,1)+(al-at).*( flux_cell_brine(:,2)).^2./vabs(:,1)+tau.*dm;%yy
Kb(:,5)=(al-at).* flux_cell_brine(:,2).* flux_cell_brine(:,3)./vabs(:,1);%yz
Kb(:,6)=at.*vabs(:,1)+(al-at).*( flux_cell_brine(:,3)).^2./vabs(:,1)+tau.*dm;%zz
Kb(isnan(Kb))=0;
rock.perm=Kb;
% water phase: half mass transmissivity due to dispersion
kb=computeTrans(G, rock);
kb(isnan(kb))=0;

%% half mass transmissivity due to dispersion in CO2-rich phase
al=state.al(2); %longitudinal dispersivity
at=state.at(2); %lateral dispersivity
dm=state.dm(2);
taubc=rock.poro.^(4/3).*(state.s0(:,2)).^(10/3);
tau=taubc;
% dispersion tensor of the CO2-rich phase
Kc(:,1)=al.*vabs(:,2)+(al-at).*( flux_cell_gas(:,1)).^2./vabs(:,2)+tau.*dm;%xx
Kc(:,2)=(al-at).* flux_cell_gas(:,1).* flux_cell_gas(:,2)./vabs(:,2);%xy
Kc(:,3)=(al-at).* flux_cell_gas(:,1).* flux_cell_gas(:,3)./vabs(:,2);%xz
Kc(:,4)=al.*vabs(:,2)+(al-at).*( flux_cell_gas(:,2)).^2./vabs(:,2)+tau.*dm;%yy
Kc(:,5)=(al-at).* flux_cell_gas(:,2).* flux_cell_gas(:,3)./vabs(:,2);%yz
Kc(:,6)=at.*vabs(:,2)+(al-at).*( flux_cell_gas(:,3)).^2./vabs(:,2)+tau.*dm;%zz
Kc(isnan(Kc))=0;
rock.perm=Kc;
% CO2-rich phase: half mass transmissivity due to dispersion
kc=computeTrans(G, rock);
kc(isnan(kc))=0;

%% half heat transmissivity

al=state.alT(1); %longitudinal dispersivity
at=state.atT(1); %lateral dispersivity
dm=state.dmT(1);
taubc=state.poro0.^(4/3).*(state.s0(:,1)).^(10/3);
tau=taubc;

% heat conductivity tensor of the water phase
Kb=zeros(size(Kb));
Kb(:,1)=al.*vabs(:,1)+(al-at).*( flux_cell_brine(:,1)).^2./vabs(:,1)+tau.*dm;%xx
Kb(:,2)=(al-at).* flux_cell_brine(:,1).* flux_cell_brine(:,2)./vabs(:,1);%xy
Kb(:,3)=(al-at).* flux_cell_brine(:,1).* flux_cell_brine(:,3)./vabs(:,1);%xz
Kb(:,4)=al.*vabs(:,1)+(al-at).*( flux_cell_brine(:,2)).^2./vabs(:,1)+tau.*dm;%yy
Kb(:,5)=(al-at).* flux_cell_brine(:,2).* flux_cell_brine(:,3)./vabs(:,1);%yz
Kb(:,6)=at.*vabs(:,1)+(al-at).*( flux_cell_brine(:,3)).^2./vabs(:,1)+tau.*dm;%zz
Kb=Kb.*state.poro0.*state.s0(:,1);
Kb(isnan(Kb))=0;
rock.perm=Kb;
% half heat transmissivity of water phase
kbT=computeTrans(G, rock);
kbT(isnan(kbT))=0;


al=state.alT(2); %longitudinal dispersivity
at=state.atT(2); %lateral dispersivity
dm=state.dmT(2);
taubc=state.poro0.^(4/3).*(state.s0(:,2)).^(10/3);
tau=taubc;
% heat conductivity tensor of the CO2-rich phase
Kc=zeros(size(Kc));
Kc(:,1)=al.*vabs(:,2)+(al-at).*( flux_cell_gas(:,1)).^2./vabs(:,2)+tau.*dm;%xx
Kc(:,2)=(al-at).* flux_cell_gas(:,1).* flux_cell_gas(:,2)./vabs(:,2);%xy
Kc(:,3)=(al-at).* flux_cell_gas(:,1).* flux_cell_gas(:,3)./vabs(:,2);%xz
Kc(:,4)=al.*vabs(:,2)+(al-at).*( flux_cell_gas(:,2)).^2./vabs(:,2)+tau.*dm;%yy
Kc(:,5)=(al-at).* flux_cell_gas(:,2).* flux_cell_gas(:,3)./vabs(:,2);%yz
Kc(:,6)=at.*vabs(:,2)+(al-at).*( flux_cell_gas(:,3)).^2./vabs(:,2)+tau.*dm;%zz
Kc=Kc.*state.poro0.*state.s0(:,2);
Kc(isnan(Kc))=0;
rock.perm=Kc;
% half heat transmissivity of CO2-rich phase
kcT=computeTrans(G, rock);
kcT(isnan(kcT))=0;


% heat conductivity of rock; n
perm_rock=state.dmT(3).*rock.poro;%(1-state.poro0).*(1-state.ss0(:,1));
% heat conductivity of hydrate;
perm_hydrate=state.dmT(4).*(rock.poro-state.poro0);
% heat conductivity of both hydrate and rock
rock.perm=perm_rock+perm_hydrate;
% half half heat transmissivity due to solid hydrate and rock
KHR=computeTrans(G, rock);
% total half heat transmissivity
KT=kbT+kcT+KHR;
end


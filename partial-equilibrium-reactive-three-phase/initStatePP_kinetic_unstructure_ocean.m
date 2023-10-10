function state = initStatePP_kinetic_unstructure_ocean(G, W, p0,rhoref,Tk,fluid,m_NaCl,rock,D_liquid,...
    D_gas,DT_liquid,DT_gas,DT_R,DT_HYD,ss0,rhomu,disco2,heat_capacity,mutj,varargin)

%
% SYNOPSIS:
%  state = initStatePP_kinetic_unstructure_ocean(G, W, p0,rhoref,Tk,fluid,m_NaCl,rock,D_liquid,...
%  D_gas,DT_liquid,DT_gas,DT_R,DT_HYD,ss0,rhomu,disco2)
%
% DESCRIPTION:
%   This function serves to initialize the system.


%
% REQUIRED PARAMETERS:
%  G  - grid system
%  W  -well information
%  p0  - initial pressure at the top of the aquifer
%  rhoref -reference density[rho_liquid,rho_gas_rho,rho_hydrate,rho_halite,rho_rock]
%  Tk  - initial temperature at the top of the aquifer
%  fluid  - fluid properties: relative permeability, retention curve, equilibrium hydration pressure 
%  m_NaCl - initial molality of NaCl in water
%  rock  - an object including the initial porosity and intrinsic permeability 
%  D_liquid- parameters for dispersion tensor of liquid phase
%  D_gas- parameters for dispersion tensor of gas phase
%  DT_liquid- parameters for heat conduction tensor of liquid phase
%  DT_gas- parameters for heat conduction tensor of gas phase
%  DT_R   -heat conductivity of rock
%  DT_HYD -heat conductivity of hydrate
%  ss0    -initial mass fractions of hydrate and halite iin the solid phase
%  rhomu  -the density and viscosity function for the fluid phases. if the  flag=0, we use the sophisticated model given by ''Numerical modeling of geological carbon  
% sequestration: enhanced dissolution in randomly heterogeneous media (https://upcommons.upc.edu/handle/2117/376077)''; if the flag=1, we use the simple linear model. 
% disco2 - an indicator for whether we consider the dissolution of CO2. if disco2=0, we do not consider the co2 dissolution,  but the value of aqueous CO2 is very small (not exactly zero)




% RETURNS:
%   state   - an object containing the initialized variables. 
%

% SEE ALSO
% D_RHOMU_kinetic_h  rhomu_p_frac_kinetic_h

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

if numel(p0) == 2
    state.pressure = repmat(p0, [G.cells.num, 1]);
end
if numel(Tk) == 1
    state.Tk = repmat(Tk, [G.cells.num, 1]);
end


if numel(rhoref) == 5
    state.rho = repmat(rhoref, [G.cells.num, 1]);
    state.rhoref=rhoref;
end
state.al=[D_liquid(1) D_gas(1)];
state.at=[D_liquid(2) D_gas(2)]; state.dm=[D_liquid(3) D_gas(3)];
state.alT=[DT_liquid(1) DT_gas(1)];
state.atT=[DT_liquid(2) DT_gas(2)]; state.dmT=[DT_liquid(3) DT_gas(3) DT_R DT_HYD];
state.heat_capacity=heat_capacity;
state.mutj=mutj;
m_NaCl=m_NaCl+linspace(0,1,G.cells.num)'.*m_NaCl.*1E-3;
s0=fluid.S(state);
state.s=s0;
state.s0=s0;
nu=[-6,-1,0,0,0,1,0;0,1,0,0,0,0,-1;0,0,-1,1,1,0,0];
omega=[1,0,0,0,0,6,0;0,1,0,0,0,1,1;0,0,1,0,1,0,0;0,0,0,1,-1,0,0];
omega_t=[1,0,0,0,0,0.710526315789474,0;0,1,0,0,0,0.289473684210526,1;0,0,1,0,1.64788732394366,0,0;0,0,0,1,-0.647887323943662,0,0];
state.nu=nu;
state.omega=omega;
state.omega_t=omega_t;
species=zeros(G.cells.num,7);
component=zeros(G.cells.num,4);
frac=zeros(G.cells.num,4);


species(:,1)=G.cells.volumes.*rock.poro.*rhoref(1).*state.s0(:,1);
species(:,4)=species(:,1).*m_NaCl.*0.023;
species(:,5)=species(:,1).*m_NaCl.*0.0355;
species(:,2)=G.cells.volumes.*rock.poro.*rhoref(2).*state.s0(:,2);%we assume the CO2 has a large quantity
species(:,3)=G.cells.volumes.*(1-rock.poro).*rhoref(4).*ss0(:,2);
species(:,6)=G.cells.volumes.*(1-rock.poro).*rhoref(3).*ss0(:,1);


state.frac=frac;
state.frac0=frac;
state.m_NaCl=ones(G.cells.num,1).*m_NaCl;
state.m_NaCl0=ones(G.cells.num,1).*m_NaCl;
state.species=species;
state.species0=species;
state.component=component;
state.component0=component;
state.poro=rock.poro;
state.poro0=rock.poro;
state.salinity=zeros(G.cells.num,1);
state.Tk0=Tk+G.cells.centroids(:,3)*0.03;
state.Tk=state.Tk0;
state.peq=fluid.peq(state);

state.rhomu=rhomu;
state.disco2=disco2;

%------------------if we have both gas and water in the beginning-----------------%

if (p0(:,2)>p0(:,1)) &&p0(:,2)<state.peq(1) &&sum(ss0)==0
    % no hydrate phase; no halite; gas and brine are at equilibrium status
    state.ss=zeros(G.cells.num,2);
    state.ss0=zeros(G.cells.num,2);
    state.rock_mass=(1-rock.poro).*G.cells.volumes.*rhoref(5).*(1-sum(state.ss,2));
    state.rock_volume=(1-rock.poro).*G.cells.volumes.*(1-sum(state.ss,2));
    state.pressure0=state.pressure;
    state.perm0=rock.perm;
    state=species_equilibrium_h_2p(state,G,rock,fluid,disco2);
    [state.rho0,state.mu0]=rhomu_p_frac_kinetic_h(state);
    [state.rho,state.mu]=deal(state.rho0,state.mu0);
    state.s0=fluid.S(state);
    state.s= state.s0;
    state.species(:,1)=G.cells.volumes.*rock.poro.*state.rho0(1).*state.s(:,1);
    state.species(:,2)=G.cells.volumes.*rock.poro.*state.rho0(2).*state.s(:,2);%we assume the CO2 has a large quantity
    state.species(:,3)=G.cells.volumes.*(1-rock.poro).*state.rho0(3).*state.ss(:,2);
    rhoeff0=effective_rho(state.rho0,state.mu0,state.s0,fluid);
    g=gravity;
    state.pressure(:,1)=state.pressure(:,1)+rhoeff0(:,1).*(G.cells.centroids(:,3)-min(G.cells.centroids(:,3))).*g(3);
    state.pressure(:,2)=state.pressure(:,2)+rhoeff0(:,1).*(G.cells.centroids(:,3)-min(G.cells.centroids(:,3))).*g(3);
    
    if ~isempty(varargin)
        if varargin{1}.type=='linear'
            dp=(varargin{1}.val(1)-varargin{1}.val(2));
            state.pressure(:,1)=state.pressure(:,1)-dp.*(G.cells.centroids(:,1)-min(G.cells.centroids(:,1)))./(max(G.cells.centroids(:,1))-min(G.cells.centroids(:,1)));
            state.pressure(:,2)=state.pressure(:,2)-dp.*(G.cells.centroids(:,1)-min(G.cells.centroids(:,1)))./(max(G.cells.centroids(:,1))-min(G.cells.centroids(:,1)));
            
        end
    end
    state.pressure0=state.pressure;
    
    
    
    
    
    
    
    for refine=1:1
        state=species_equilibrium_h_2p(state,G,rock,fluid,disco2);
        
        
        state.s0=fluid.S(state);
        state.s=state.s0;
        [state.rho0,state.mu0]=rhomu_p_frac_kinetic_h(state);
        state.rho=state.rho0;
        state.mu=state.mu0;
        state.species0=state.species;
        state.component0=state.component;
        state.frac0=state.frac;
        state.salinity0=state.salinity;
        state.ss0=state.ss;
        state=species_equilibrium_h_2p(state,G,rock,fluid,disco2);
        
        %%
        state.poro0=state.poro;
        state.ss0=state.ss;
        state.species(:,3)=state.ss(:,2).*G.cells.volumes.*(1-state.poro).*state.rho(:,3);
        state.species(:,2)=state.rho(:,2).*state.s(:,2).*rock.poro(:).*G.cells.volumes(:);
        state.species(:,6)=state.rho(:,3).*(1-rock.poro(:)).*G.cells.volumes(:).*state.ss(:,1);
        state.species(:,[1,4:5,7])=state.rho(:,1).*state.s(:,1).*rock.poro(:).*G.cells.volumes(:).*state.species(:,[1,4:5,7])./(sum(state.species(:,[4:5,7]),2)+state.species(:,1));
        for i= 1:G.cells.num
            speciesw=[state.species(i,1);0;0;0];
            state.frac(i,:)=(state.omega_t(:,[4:5,7])*state.species(i,[4:5,7])'+speciesw)'./(sum(state.species(i,:))-state.species(i,2)-state.species(i,3)-state.species(i,6));
            state.salinity(i,1)=((sum(state.species(i,[4:5,7]),2))./(sum(state.species(i,[4:5,7]),2)+state.species(i,1)))';
            state.m_NaCl(i)=0.5.*(state.species(i,4)/0.023+state.species(i,5)./0.0355)./state.species(i,1);
            state.m_NaCl0(i)=state.m_NaCl(i);
        end
        
        state.species0=state.species;
        state.frac0=state.frac;
        state.salinity0=state.salinity;
        [d_rho,~]=D_RHOMU_kinetic_h(state);
        [state.h0,state.u0,~,~]=heat_h(state,d_rho);
        state.h=state.h0;
        state.u=state.u0;
        
        for i=1:G.cells.num
            
            state.component(i,:)=state.omega_t*state.species(i,:)';
            
        end
        
    end
end

%------------------if we have only water in the beginning-----------------%

if (p0(:,2)<=p0(:,1)) % one phase
    
    state.ss=zeros(G.cells.num,2);
    state.ss0=zeros(G.cells.num,2);
    state.rock_mass=(1-rock.poro).*G.cells.volumes.*rhoref(5).*(1-sum(state.ss,2));
    state.rock_volume=(1-rock.poro).*G.cells.volumes.*(1-sum(state.ss,2));
    state.pressure0=state.pressure;
    state.perm0=rock.perm;
    state=species_equilibrium_h_2p(state,G,rock,fluid,disco2);
    [state.rho0,state.mu0]=rhomu_p_frac_kinetic_h(state);
    [state.rho,state.mu]=deal(state.rho0,state.mu0);
    state.s0=fluid.S(state);
    state.s= state.s0;
    state.species(:,1)=G.cells.volumes.*rock.poro.*state.rho0(1).*state.s(:,1);
    state.species(:,2)=G.cells.volumes.*rock.poro.*state.rho0(2).*state.s(:,2);%we assume the CO2 has a large quantity
    state.species(:,3)=G.cells.volumes.*(1-rock.poro).*state.rho0(3).*state.ss(:,2);
    rhoeff0=effective_rho(state.rho0,state.mu0,state.s0,fluid);
    g=gravity;
    state.pressure(:,1)=state.pressure(:,1)+rhoeff0(:,1).*(G.cells.centroids(:,3)-min(G.cells.centroids(:,3))).*g(3);
    state.pressure(:,2)=state.pressure(:,2);
    state.pressure0=state.pressure;
     
    
    
    for refine=1:1
        state=species_equilibrium_h_2p(state,G,rock,fluid,disco2);      
        state.s0=fluid.S(state);
        state.s=state.s0;
        [state.rho0,state.mu0]=rhomu_p_frac_kinetic_h(state);
        state.rho=state.rho0;
        state.mu=state.mu0;
        state.species0=state.species;
        state.component0=state.component;
        state.frac0=state.frac;
        state.salinity0=state.salinity;
        state.ss0=state.ss;
        state=species_equilibrium_h_2p(state,G,rock,fluid,disco2);
        
        %%
        state.poro0=state.poro;
        state.ss0=state.ss;
        state.species(:,3)=state.ss(:,2).*G.cells.volumes.*(1-state.poro).*state.rho(:,3);
        state.species(:,2)=state.rho(:,2).*state.s(:,2).*rock.poro(:).*G.cells.volumes(:);
        state.species(:,6)=state.rho(:,3).*(1-rock.poro(:)).*G.cells.volumes(:).*state.ss(:,1);
        state.species(:,[1,4:5,7])=state.rho(:,1).*state.s(:,1).*rock.poro(:).*G.cells.volumes(:).*state.species(:,[1,4:5,7])./(sum(state.species(:,[4:5,7]),2)+state.species(:,1));
        for i= 1:G.cells.num
            speciesw=[state.species(i,1);0;0;0];
            state.frac(i,:)=(state.omega_t(:,[4:5,7])*state.species(i,[4:5,7])'+speciesw)'./(sum(state.species(i,:))-state.species(i,2)-state.species(i,3)-state.species(i,6));
            state.salinity(i,1)=((sum(state.species(i,[4:5,7]),2))./(sum(state.species(i,[4:5,7]),2)+state.species(i,1)))';
            state.m_NaCl(i)=0.5.*(state.species(i,4)/0.023+state.species(i,5)./0.0355)./state.species(i,1);
            state.m_NaCl0(i)=state.m_NaCl(i);
        end
        
        state.species0=state.species;
        state.frac0=state.frac;
        state.salinity0=state.salinity;
        [d_rho,~]=D_RHOMU_kinetic_h(state);
        [state.h0,state.u0,~,~]=heat_h(state,d_rho);
        state.h=state.h0;
        state.u=state.u0;
        
        for i=1:G.cells.num
            
            state.component(i,:)=state.omega_t*state.species(i,:)';
            
        end
        
    end
    
end



if ~isempty(W)
    state.wellSol = initWellSol(W, p0);
end
state.species_name={'H2O','CO2(g)','NaCl(s)','Na+','Cl-','HYD(s)','CO2(aq)'};
end

function rhoeff = effective_rho(rho0,mu0,s0, fluid)


kr        = fluid.relperm(s0);
mob    = bsxfun(@rdivide, kr, mu0);
totmob = sum(mob, 2);
rhoeff = sum(bsxfun(@times, mob, rho0(:,1:2)), 2) ./ totmob;
end



function state=species_equilibrium_h_2p(state,G,rock,fluid,disco2)


species_t=zeros(G.cells.num,7);
for i=1:G.cells.num
    species_t(i,:)=reaction_equilibrium_h_2p(state.pressure(i,1),state.pressure(i,2),state.Tk(i),state.species(i,:),disco2);
end

for j=1:1
    [state.rho,~]=rhomu_p_frac_kinetic_h(state);
    
    dt=1e2;
    
    for i=1:G.cells.num
        % component_t=state.omega_t*state.species(i,:)';
        
        state.species(i,:)=species_t(i,:);
        state.species(i,2)=state.rho(i,2).*state.s(i,2).*rock.poro(i).*G.cells.volumes(i);
        state.species(i,3)=state.rho(i,3).*(1-rock.poro(i)).*G.cells.volumes(i).*state.ss(i,2);
        state.species(i,6)=state.rho(i,3).*(1-rock.poro(i)).*G.cells.volumes(i).*state.ss(i,1);
        state.species(i,1)=state.rho(i,1).*state.s(i,1).*rock.poro(i).*G.cells.volumes(i).*state.species(i,1)./(sum(state.species(i,4:5),2)+state.species(i,1)+state.species(i,7));
        state.species(i,4)=state.species(i,1).*state.m_NaCl(i).*0.023;
        state.species(i,5)=state.species(i,1).*state.m_NaCl(i).*0.0355;
        state.species(i,[1,4:5,7])=state.rho(i,1).*state.s(i,1).*rock.poro(i).*G.cells.volumes(i).*state.species(i,[1,4:5,7])./(sum(state.species(i,[4:5,7]),2)+state.species(i,1));
        
        state.m_NaCl(i)=0.5.*(state.species(i,4)/0.023+state.species(i,5)./0.0355)./state.species(i,1);
        state.m_NaCl0(i)=state.m_NaCl(i);       
        state.component(i,:)=state.omega_t*state.species(i,:)';
        
    end
    state.component0=state.component;
    state.species0=state.species;
    state=kinetic_reaction_hc0(state,G.cells.volumes,fluid,dt);
    error= (state.component-state.component0)./(abs(state.component0)+abs(state.component));
    max(error);
    mean(error);
    for  i=1:G.cells.num
        speciesw=[state.species(i,1);0;0;0];
        state.frac(i,:)=(state.omega_t(:,[4:5,7])*state.species(i,[4:5,7])'+speciesw)'./(sum(state.species(i,:))-state.species(i,2)-state.species(i,3)-state.species(i,6));
        state.salinity(i,1)=((sum(state.species(i,[4:5,7]),2))./(sum(state.species(i,[4:5,7]),2)+state.species(i,1)))';
    end
end



end

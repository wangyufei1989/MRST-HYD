function state = initStatePP_kinetic_unstructure_hc(G, W, p0,rhoref,Tk,fluid,m_NaCl,rock,D_liquid,...
    D_gas,DT_liquid,DT_gas,DT_R,DT_HYD,ss0)
% initialize the system.
%
% SYNOPSIS:
%  state = initStatePP_kinetic1(G, W, p0,rhoref,muref,Tk,fluid,m_NaCl,rock,property)
%
% DESCRIPTION:
%   This function serves to initialize the system.


%
% REQUIRED PARAMETERS:
%   G  - grid system
%  W  -well information
%  p0  - initial pressure at the top of the aquifer
%   s0  - initial saturation at the top of the aquifer
% frac0  - intial mass fraction at the top of the aquifer
%rhoref, muref  - reference density and viscosity; not used if the density
%and viscosity change with pressure.
%fluid  - relative permeability and retention curves;
% m_NaCl,Tk  -  (constant) NaCl concentration and (constant) temperature 
%rock  - an object including the initial porosity and intrinsic permeability 
%property  - indicates what phases and reactions do we have.


% RETURNS:
%   state   - an object containing the initialized variables. 
%

%{

This file is part of mrst_co2 based on MRST.

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

%
property=5;
state.Tk0=Tk+G.cells.centroids(:,3)*0.035;
%state.Tk0=Tk+G.cells.centroids(:,3)*0;
state.Tk=state.Tk0;
state.peq=fluid.peq(state);



 if (p0(:,2)>p0(:,1)) &&p0(:,2)<state.peq(1) && property==5 &&sum(ss0)==0% two phases;
     % no hydrate phase; no halite; gas and brine are at equilibrium status  
     state.ss=zeros(G.cells.num,2);
     state.ss0=zeros(G.cells.num,2);
     state.rock_mass=(1-rock.poro).*G.cells.volumes.*rhoref(5).*(1-sum(state.ss,2));
state.rock_volume=(1-rock.poro).*G.cells.volumes.*(1-sum(state.ss,2));
     state.pressure0=state.pressure;
     state.perm0=rock.perm;
       state=species_equilibrium_h_2p(state,G,rock,fluid);
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
       eric=0;
       
       if eric==1
            state.pressure(:,1)=state.pressure(:,1)-2e6.*(G.cells.centroids(:,1)-min(G.cells.centroids(:,1)))./(max(G.cells.centroids(:,1))-min(G.cells.centroids(:,1)));
       state.pressure(:,2)=state.pressure(:,2)-2e6.*(G.cells.centroids(:,1)-min(G.cells.centroids(:,1)))./(max(G.cells.centroids(:,1))-min(G.cells.centroids(:,1)));   
       end
      state.pressure0=state.pressure;

      
       

      
      
      
      for refine=1:1
      state=species_equilibrium_h_2p(state,G,rock,fluid);
 

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
      state=species_equilibrium_h_2p(state,G,rock,fluid);
    
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
  [d_rho,d_mu]=D_RHOMU_kinetic_h(state);    
 [state.h0,state.u0,~,~]=heat_h(state,d_rho);
   state.h=state.h0;
   state.u=state.u0;
   
      for i=1:G.cells.num

state.component(i,:)=state.omega_t*state.species(i,:)';

      end
   
 end
      end
       
  %%     
       
        if (p0(:,2)<=p0(:,1))  && property==5% one phase
  
      [rho0,mu0]=rhomu_p_frac_kinetic(state);
       s0=fluid.S(state);
       rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        g=gravity;
       state.pressure(:,1)=state.pressure(:,1)+rhoeff0(:,1).*(G.cells.centroids(:,3)-min(G.cells.centroids(:,3))).*g(3);
           state.pressure(:,2)=state.pressure(:,2)+rhoeff0(:,1).*(G.cells.centroids(:,3)-min(G.cells.centroids(:,3))).*g(3);   
     
          [rho0,mu0]=rhomu_p_frac_kinetic(state);
          [state.rho,state.mu]=rhomu_p_frac_kinetic(state);
               state=species_equilibrium(state,G,rock);
           s0=fluid.S(state);
       rhoeff0=effective_rho(rho0,mu0,s0,fluid);
       

for refine=1:50
       %  state=species_equilibrium(state,G,rock);
end

 
        state.s0=fluid.S(state);
        state.s=state.s0;
      [state.rho0,state.mu0]=rhomu_p_frac_kinetic(state);
      state.rho=state.rho0;
      state.mu=state.mu0;
      state.species0=state.species;
      state.component0=state.component;
      state.frac0=state.frac;
      state.salinity0=state.salinity;
    
  end
 


%state=species_equilibrium(state,G,rock);

   if ~isempty(W)
      state.wellSol = initWellSol(W, p0);
   end
end

function rhoeff = effective_rho(rho0,mu0,s0, fluid)
   
 
   kr        = fluid.relperm(s0);

   mob    = bsxfun(@rdivide, kr, mu0);
   totmob = sum(mob, 2);
   rhoeff = sum(bsxfun(@times, mob, rho0(:,1:2)), 2) ./ totmob;
end

function state=species_equilibrium(state,G,rock)


for i=1:G.cells.num
   % component_t=state.omega_t*state.species(i,:)';
species_t=reaction_equilibrium(state.pressure(i,1),state.pressure(i,2),state.nu,state.omega,state.species(i,:));
state.species(i,:)=species_t';
state.species(i,2)=state.rho(i,2).*state.s(i,2).*rock.poro(i).*G.cells.volumes(i);
state.species(i,3)=state.rho(i,3).*(1-rock.poro(i)).*G.cells.volumes(i);
state.species(i,1)=state.rho(i,1).*state.s(i,1).*rock.poro(i).*G.cells.volumes(i).*state.species(i,1)./(sum(state.species(i,4:end),2)+state.species(i,1)*state.m_NaCl*0.05844+state.species(i,1));
state.species(i,4:end)=state.rho(i,1).*state.s(i,1).*rock.poro(i).*G.cells.volumes(i).*state.species(i,4:end)./(sum(state.species(i,4:end),2)+state.species(i,1)*state.m_NaCl*0.05844+state.species(i,1));

state.component(i,:)=state.omega_t*species_t;
speciesw=[state.species(i,1);0;0;0];
state.frac(i,:)=(state.omega_t(:,4:end)*state.species(i,4:end)'+speciesw)'./(sum(state.species(i,:))-state.species(i,2)-state.species(i,3)+state.species(i,1)*state.m_NaCl*0.05844);
state.salinity(i)=(sum(state.species(i,4:end),2)+state.species(i,1)*state.m_NaCl*0.05844)./(sum(state.species(i,4:end),2)+state.species(i,1)*state.m_NaCl*0.05844+state.species(i,1));
end

     

end

function state=species_equilibrium_h_2p(state,G,rock,fluid)


species_t=zeros(G.cells.num,7);
for i=1:G.cells.num
 species_t(i,:)=reaction_equilibrium_h_2p(state.pressure(i,1),state.pressure(i,2),state.Tk(i),state.species(i,:));   
end

for j=1:1
   [state.rho,~]=rhomu_p_frac_kinetic_h(state);

dt=1e6;
  
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
%state.species(i,4)=state.species(i,1).*state.m_NaCl(i).*0.023;
%state.species(i,5)=state.species(i,1).*state.m_NaCl(i).*0.0355;

state.component(i,:)=state.omega_t*state.species(i,:)';

   end
   state.component0=state.component;
state.species0=state.species;

 state=kinetic_reaction_hc(state,G.cells.volumes,fluid,dt);
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


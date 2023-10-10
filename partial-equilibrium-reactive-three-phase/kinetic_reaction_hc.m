function state=kinetic_reaction_hc(state,V,fluid,DT,rock)

% speciation calculation based on the kinetic reaction given in section called
%Partial-Equilibrium Reactive Multi-Component Three-Phase Flow Mode in the user'guides
%
% SYNOPSIS:
% state=kinetic_reaction(state,V,DT)
%
% DESCRIPTION:
%   This function serves to calculate the mass fraction of each species and components. 
%Description is already included in the user's guide. 

%
% REQUIRED PARAMETERS:
% state   -state varibles
% V         - cell volume
% DT      -time step 

% RETURNS:
%   state   - contains the mass fractions of all species and components at equilibrium state; 
%                and also their derivatives over gas pressure. 

% SEE ALSO:
%  kr3p

%{

This file is part of mrst_co2h based on MRST.

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

%% main file for kinetic reaction

state.species(:,6)=state.species0(:,6);% the hydrate mass is equal to that in the previous time step to avoid repeated calculation
speciestem1=state.species;% this is used for the calculation of reaction when pressure is reduced
speciestem2=state.species;% this is used for the calculation of reaction when the temperature is increased
component=state.component;% component does not change during the reaction
frac=state.frac;% initialize the mass fraction of each component in water; the value can be any (e.g., zero); just for generate variable space.
salinity=zeros(size(state.pressure,1));% initialize the salinity of water; the value can be any; just for generate variable space.
mc0=zeros(size(state.pressure,1));% initialize the mass fraction of CO2 in water; the value can be any; just for generate variable space.
mc=zeros(size(state.pressure,1));% initialize the mass fraction of CO2 in water; the value can be any; just for generate variable space.
state.dsp=zeros(size(V));% generate the variable space for d_salinity/d_pg
state.dcp=zeros(size(V));% generate the variable space for d_X_l^C/d_pg; X_l^C is mass fraction of CO2 in water
  dpp=-1e0;% pressure perturbation for calculation d (...)/dp
  dTT=1E-2;% temperature perturbation for calculation d (...)/dp
 disco2=state.disco2;
for i=1:size(state.pressure,1)
    if i==781
     % pause
    end
    
    
k0=state.kinetic_rate;
Sr=V(i).*state.poro0(i).^1.5.*(2.*state.permr(i))^(-.5).*k0;% k=Sr*exp(-Ea/RT)*(p-p_eq); 
%----------------------------------------------------------------------------------------
[ nt,comp,kr,dr_k_pg]=kinetic_reaction_hc_sub(state.pressure(i,1),state.pressure(i,2),state.peq(i),state.nu,state.omega,state.omega_t,state.species(i,:),state.component(i,:),state.Tk(i),Sr,DT,state.rock_mass(i),...
         state.pressure(i,2),state.Ea1,state.Ea2,state.nr_chemical(i),disco2);% speciation; partial equilibrium reaction of the chemical system

%---------------------------------------------------------------------------------------- 
 state.dr_k_pg(i)=dr_k_pg;
 state.species(i,:)=nt;% masses of all species in the i-th cell  
 state.component(i,:)=comp';% actually, this is not necessary, because the component does not change during reaction; this is only for checking if the chemical system is correct; 
 state.m_NaCl(i)=0.5.*(nt(4)./0.023+nt(5)./0.0355)./nt(1); %the molality of NaCl in water
 state.ss(i,2)=state.species(i,3)./(state.species(i,3)+state.species(i,6)+state.rock_mass(i));% mass fraction of halite in solid phase; currently, this is not used
 state.ss(i,1)=state.species(i,6)./(state.species(i,3)+state.species(i,6)+state.rock_mass(i));% mass fraction of hydrate in solid phase; currently, this is not used
 state.r_k(i,1)=kr;% kinetic rate of hydrate formation (+) dissociation (-)
 speciesw=[state.species(i,1);0;0;0];% water mass
 state.frac(i,:)=(state.omega_t(:,[4:5,7])*state.species(i,[4:5,7])'+speciesw)'./(sum(state.species(i,:))-state.species(i,2)-state.species(i,3)-state.species(i,6));% mass fraction of each component in water
 state.poroh(i)=state.species(i,6)./state.rhoref(3)/V(i);%volume fraction of hydrate; phi_hydrate
 state.porona(i,1)=state.species(i,3)./state.rhoref(4)/V(i);%volume fraction of halite; phi_halite; negligible
 state.salinity(i)=(sum(state.species(i,[4:5,7]),2))./(sum(state.species(i,[4:5,7]),2)+state.species(i,1));% salinity
 mc0(i)=state.species(i,7)./(sum(state.species(i,[4:5,7]),2)+state.species(i,1)); %mass fraction of CO2 in water 
    

end

%%
state.pressure(:,2)=state.pressure(:,2)+dpp; % increase the pressure by dpp to calculate d(...)/dp; here dpp can be posirive or negative, as long as it is small
component(:,2)=component(:,2)*1; % here, we increase the mass of CO2 component to ensure that the CO2 supply is enough; for instance, if the co2 mass is 10 kg and it is in aqueous state,
%but changing the gas pressure can make the dissolved mass increase; therefore, we need more co2 to assure that the co2 is abundant 

state.dsh=zeros(size(state.ss));% generate varaible space for derivative of hydrate or halite mass fraction in solid phase over gas pressure; curretnly, this is not used
species=zeros(size(state.species));% generate variable space for species
 ss=zeros(size(state.ss));% generate variable space for mass fraction of halite/hydrate in solid phase; currently, this is not used
for  i=1:size(state.pressure,1)
    k0=state.kinetic_rate;
Sr=V(i).*state.poro0(i).^1.5.*(2.*state.permr(i))^(-.5).*k0;% k=Sr*exp(-Ea/RT)*(p-p_eq); 
 %----------------------------------------------------------------------------------------
 [ nt,comp,kr,dr_k_pg]=kinetic_reaction_hc_sub(state.pressure(i,1),state.pressure(i,2),state.peq(i),state.nu,state.omega,state.omega_t,speciestem1(i,:),component(i,:),state.Tk(i),Sr,DT,state.rock_mass(i),...
        state.pressure(i,2),state.Ea1,state.Ea2,state.nr_chemical(i),disco2);% speciation; partial equilibrium reaction of the chemical system
 %----------------------------------------------------------------------------------------  
species(i,:)=nt;% mass of all species in i-th cell
specieswt=[species(i,1);0;0;0];% mass of water species
frac(i,:)=(state.omega_t(:,[4:5,7])*species(i,[4:5,7])'+specieswt)'./(sum(species(i,:))-species(i,2)-species(i,3)-species(i,6));%mass fraction of each component in water
salinity(i)=(sum(species(i,[4:5,7]),2))./(sum(species(i,[4:5,7]),2)+species(i,1));% salinity
mc(i)=species(i,7)./(sum(species(i,[4:5,7]),2)+species(i,1));%mass fraction of CO2 in water 
state.dcp(i)=(mc(i)-mc0(i))'./dpp;% d_mc/d_p
state.dsp(i)=(salinity(i)-state.salinity(i))'./dpp;%d_salinity/dp  
state.dfrac(i,:)=(frac(i,:)-state.frac(i,:))./dpp;%derivative of mass fraction of component in water over pressure
ss(i,2)=species(i,3)./(species(i,3)+species(i,6)+state.rock_mass(i));% mass fraction of halite in solid phase; currently, this is not used
ss(i,1)=species(i,6)./(species(i,3)+species(i,6)+state.rock_mass(i));% mass fraction of hydrate in solid phase; currently, this is not used
state.dsh(i,1)=(ss(i,1)-state.ss(i,1))./dpp';% derivatie of mass fraction of hydrate in solid phase over pressure; currently, this is not used
state.dsh(i,2)=(ss(i,2)-state.ss(i,2))./dpp';% derivatie of mass fraction of halite in solid phase over pressure; currently, this is not used
state.dphip(i,1)=-(((species(i,6)-state.species(i,6))./state.rhoref(3)+(species(i,3)-state.species(i,3))./state.rhoref(4))./V(i)./dpp)';% derivative of porosity over pressure
state.dphihp(i,1)=-(((species(i,6)-state.species(i,6))./state.rhoref(3))./V(i)./dpp)';% derivative of porosity over pressure (neglect the halite)
state.dspecies(i,:)=(species(i,:)-state.species(i,:))./dpp;% derivatives of species over pressure; warning CO2 gas species cannot be used
    
end
%%
state.pressure(:,2)=state.pressure(:,2)-dpp;% the gas pressure is reovered
state.Tk=state.Tk+dTT;% perturb the temperature
state.peq=fluid.peq(state);% perturbation of the temperature can influence the equilibrium pressure for hydration 
state.dshT=zeros(size(state.ss));%generate varaible space for derivative of hydrate or halite mass fraction in solid phase over temperature; curretnly, this is not used
for  i=1:size(state.pressure,1)
 k0=state.kinetic_rate;
Sr=V(i).*state.poro0(i).^1.5.*(2.*state.permr(i))^(-.5).*k0;% k=Sr*exp(-Ea/RT)*(p-p_eq);    
   %----------------------------------------------------------------------------------------  

     [ nt,comp,kr,dr_k_pg]=kinetic_reaction_hc_sub(state.pressure(i,1),state.pressure(i,2),state.peq(i),state.nu,state.omega,state.omega_t,speciestem2(i,:),component(i,:),state.Tk(i),Sr,DT,state.rock_mass(i),...
        state.pressure(i,2),state.Ea1,state.Ea2,state.nr_chemical(i),disco2);% speciation; partial equilibrium reaction of the chemical system
  %----------------------------------------------------------------------------------------  
species(i,:)=nt;% mass of all species in i-th cell
specieswt=[species(i,1);0;0;0];% mass of water species
frac(i,:)=(state.omega_t(:,[4:5,7])*species(i,[4:5,7])'+specieswt)'./(sum(species(i,:))-species(i,2)-species(i,3)-species(i,6));%mass fraction of each component in water
salinity(i)=(sum(species(i,[4:5,7]),2))./(sum(species(i,[4:5,7]),2)+species(i,1));% salinity
mc(i)=species(i,7)./(sum(species(i,[4:5,7]),2)+species(i,1));%mass fraction of CO2 in water 
state.dcpT(i)=(mc(i)-mc0(i))'./dTT;%d_mc/dT
state.dspT(i)=(salinity(i)-state.salinity(i))'./dTT;%d_salinity/dT
state.dfracT(i,:)=(frac(i,:)-state.frac(i,:))./dTT;% derivative of mass fraction of component in water over temperature
ss(i,2)=species(i,3)./(species(i,3)+species(i,6)+state.rock_mass(i));% mass fraction of halite in solid phase; currently, this is not used
ss(i,1)=species(i,6)./(species(i,3)+species(i,6)+state.rock_mass(i));% mass fraction of hydrate in solid phase; currently, this is not used      
state.dshT(i,1)=(ss(i,1)-state.ss(i,1))./dTT'; % derivatie of mass fraction of hydrate in solid phase over temperature; currently, this is not used
state.dshT(i,2)=(ss(i,2)-state.ss(i,2))./dTT'; % derivatie of mass fraction of halite in solid phase over temperature; currently, this is not used
state.dphiT(i,1)=-(((species(i,6)-state.species(i,6))./state.rhoref(3)+(species(i,3)-state.species(i,3))./state.rhoref(4))./V(i)./dTT)';% derivative of porosity over temperature
state.dphihT(i,1)=-(((species(i,6)-state.species(i,6))./state.rhoref(3))./V(i)./dTT)';% derivative of porosity over temperature (neglect the halite)
state.dspeciesT(i,:)=(species(i,:)-state.species(i,:))./dTT;% derivatives of species over pressure; warning CO2 gas species cannot be used
    
end
state.Tk=state.Tk-dTT;% recover T

end
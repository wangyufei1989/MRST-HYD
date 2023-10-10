function state=kinetic_reaction_hc0(state,V,fluid,DT)

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

%% main file for kinetic reaction

state.species(:,6)=state.species0(:,6);
speciestem1=(state.species);
speciestem2=state.species;
component=state.component;
frac=state.frac;
salinity=zeros(size(state.pressure,1));
state.dsp=zeros(size(V));
state.dcp=zeros(size(V));
  dpp=-1e2;
  dTT=1E-2;
 
for i=1:size(state.pressure,1)
k0=1e-14*1e3;

if state.pressure(i,2)<state.peq(i)
      Sr=V(i).*state.poro0(i).^1.5.*(2.*state.perm0(i)).^(-.5).*k0;%k0*Gamma_r*Ah; where Ah=phi^1.5/(2kappa)^.5; Gamma_r=hydrate saturation 
else
     Sr=V(i).*state.poro0(i).^1.5.*(2.*state.perm0(i))^(-.5).*k0; %k0*Gamma_r*Ah; where Ah=phi^1.5/(2kappa)^.5; Gamma_r=s_l*s_g
% k=Sr*exp(-Ea/RT)*(p_eq-p); so  we still need exp(-Ea/RT)*(p_eq-p)
end



   %p1=ones(size(state.pressure,1),1).*1e7;
   sr_chemical=1e-6;ksal=10^(1.582);
   if (state.s(i,1)<sr_chemical||state.s(i,2)<sr_chemical) && state.pressure0(i,2)>state.peq(i)  
      Sr=0;  
      ksal=10^(1.582)*10;
   end
     [ nt,comp,kr,dr_k_pg]=kinetic_reaction_h_sub(state.pressure(i,1),state.pressure(i,2),state.peq(i),state.nu,state.omega,state.omega_t,state.species(i,:),state.component(i,:),state.Tk(i),Sr,DT,state.rock_mass(i),...
         ksal,state.pressure(i,2),state.disco2);
  
   state.dr_k_pg(i)=dr_k_pg;
   state.species(i,:)=nt;
 
     state.component(i,:)=comp';
  state.m_NaCl(i)=0.5.*(nt(4)./0.023+nt(5)./0.0355)./nt(1); 
  
    state.ss(:,2)=state.species(:,3)./(state.species(:,3)+state.species(:,6)+state.rock_mass);
 state.ss(:,1)=state.species(:,6)./(state.species(:,3)+state.species(:,6)+state.rock_mass);
%    
  state.kr(i,1)=kr;
  speciesw=[state.species(i,1);0;0;0];
  state.frac(i,:)=(state.omega_t(:,[4:5,7])*state.species(i,[4:5,7])'+speciesw)'./(sum(state.species(i,:))-state.species(i,2)-state.species(i,3)-state.species(i,6));

  state.poroh(i,1)=state.species(i,6)./state.rhoref(3)/V(i);
   state.porona(i,1)=state.species(i,3)./state.rhoref(4)/V(i);
  
state.salinity(i)=(sum(state.species(i,[4:5,7]),2))./(sum(state.species(i,[4:5,7]),2)+state.species(i,1));
 mc0=state.species(i,7)./(sum(state.species(i,[4:5,7]),2)+state.species(i,1));  
    

end

%%
state.pressure(:,2)=state.pressure(:,2)+dpp;
% component(:,1)=component(:,1)+component(:,2).*.7105./.2895;
% component(:,2)=component(:,2)+component(:,2);


 state.peq=fluid.peq(state);
     state.dsh=zeros(size(state.ss));
for  i=1:size(state.pressure,1)
    

     [ nt,comp,kr,dr_k_pg]=kinetic_reaction_h_sub(state.pressure(i,1),state.pressure(i,2),state.peq(i),state.nu,state.omega,state.omega_t,speciestem1(i,:),component(i,:),state.Tk(i),Sr,DT,state.rock_mass(i),...
         ksal,state.pressure(i,2),state.disco2);
   
     species(i,:)=nt;

   component(i,:)=comp';
  
   
 
 specieswt=[species(i,1);0;0;0];
frac(i,:)=(state.omega_t(:,[4:5,7])*species(i,[4:5,7])'+specieswt)'./(sum(species(i,:))-species(i,2)-species(i,3)-species(i,6));


salinity(i)=(sum(species(i,[4:5,7]),2))./(sum(species(i,[4:5,7]),2)+species(i,1));


mc=species(i,7)./(sum(species(i,[4:5,7]),2)+species(i,1));
state.dcp(i)=(mc-mc0)'./dpp;
state.dsp(i)=(salinity(i)-state.salinity(i))'./dpp;

state.dfrac(i,:)=(frac(i,:)-state.frac(i,:))./dpp;

 ss=state.ss;
   ss(i,2)=species(i,3)./(species(i,3)+species(i,6)+state.rock_mass(i));
       ss(i,1)=species(i,6)./(species(i,3)+species(i,6)+state.rock_mass(i));
  
      state.dsh(i,1)=(ss(i,1)-state.ss(i,1))./dpp';
      state.dsh(i,2)=(ss(i,2)-state.ss(i,2))./dpp';
  

state.dphip(i,1)=-(((species(i,6)-state.species(i,6))./state.rhoref(3)+(species(i,3)-state.species(i,3))./state.rhoref(4))./V(i)./dpp)';

state.dphihp(i,1)=-(((species(i,6)-state.species(i,6))./state.rhoref(3))./V(i)./dpp)';

  state.dspecies(i,:)=(species(i,:)-state.species(i,:))./dpp;
    
end
%%
state.pressure(:,2)=state.pressure(:,2)-dpp;
state.Tk=state.Tk+dTT;

 state.peq=fluid.peq(state);
 state.dshT=zeros(size(state.ss));
for  i=1:size(state.pressure,1)
    

     [ nt,comp,kr,dr_k_pg]=kinetic_reaction_h_sub(state.pressure(i,1),state.pressure(i,2),state.peq(i),state.nu,state.omega,state.omega_t,speciestem2(i,:),component(i,:),state.Tk(i),Sr,DT,state.rock_mass(i),...
         ksal,state.pressure(i,2),state.disco2);
   
     species(i,:)=nt;

   component(i,:)=comp';
  
   
  specieswt=[species(i,1);0;0;0];
frac(i,:)=(state.omega_t(:,[4:5,7])*species(i,[4:5,7])'+specieswt)'./(sum(species(i,:))-species(i,2)-species(i,3)-species(i,6));

salinity(i)=(sum(species(i,[4:5,7]),2))./(sum(species(i,[4:5,7]),2)+species(i,1));


mc=species(i,7)./(sum(species(i,[4:5,7]),2)+species(i,1));
state.dcpT(i)=(mc-mc0)'./dTT;
state.dspT(i)=(salinity(i)-state.salinity(i))'./dTT;

state.dfracT(i,:)=(frac(i,:)-state.frac(i,:))./dTT;

 ss=state.ss;
   ss(i,2)=species(i,3)./(species(i,3)+species(i,6)+state.rock_mass(i));
       ss(i,1)=species(i,6)./(species(i,3)+species(i,6)+state.rock_mass(i));
      
      state.dshT(i,1)=(ss(i,1)-state.ss(i,1))./dTT';
      state.dshT(i,2)=(ss(i,2)-state.ss(i,2))./dTT';
  
  

state.dphiT(i,1)=-(((species(i,6)-state.species(i,6))./state.rhoref(3)+(species(i,3)-state.species(i,3))./state.rhoref(4))./V(i)./dTT)';

state.dphihT(i,1)=-(((species(i,6)-state.species(i,6))./state.rhoref(3))./V(i)./dTT)';

  state.dspeciesT(i,:)=(species(i,:)-state.species(i,:))./dTT;
    
end
rf=0;
rf1=0;
 state.dphiT=rf.*state.dphiT;
 state.dphip=rf.*state.dphip;
  state.dspecies=rf.*state.dspecies;
 state.dspeciesT=rf.*state.dspeciesT;
 state.dfrac=rf1.*state.dfrac;
  state.dfracT=rf1.*state.dfracT;
 state.dsh=rf.*state.dsh;
 state.dshT=rf.*state.dshT;
    state.dsp=rf.*state.dsp;
 state.dspT=rf.*state.dspT;
   state.dcp=rf1.*state.dcp;
 state.dcpT=rf1.*state.dcpT;
 
 
 
 
state.Tk=state.Tk-dTT;


end

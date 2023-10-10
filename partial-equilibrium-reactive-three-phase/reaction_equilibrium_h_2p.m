function n=reaction_equilibrium_h_2p(Pl,Pg,Tk,n,disco2)

%
% SYNOPSIS:
%  n=reaction_equilibrium_h_2p(Pl,Pg,Tk,n,disco2)
%
% DESCRIPTION:
%   This function serves to calculate the mass fraction of CO2 component in water. 


%
% REQUIRED PARAMETERS:
% Pl       -liquid pressure
%Pg      - gas pressure
% nu   -stoichiometric matrix
% G      -kernel matrix  based on molar abundance
% nu      - stoichiometric matrix 
% n         -mass of the species
% h         - mass of the components



% RETURNS:
%   n   -  mass of species 



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
%% change pressure from [pa] to [bar]
Pg=Pg.*1e-5;
Pl=Pl.*1e-5;
%% model for activity coefficient of aqueous CO2
gammac='duan_sun2003';
if gammac=='drummond1981'
% Drummond 1981 model; see ''Numerical modeling of geological carbon  
% sequestration: enhanced dissolution in randomly heterogeneous media (https://upcommons.upc.edu/handle/2117/376077)'' 
ap1=-1.0312;
ap2=1.2806e-3;
ap3=255.9;
ap4=0.4445;
ap5=-1.606e-3;
m_na=n(4)./0.023./n(1);
m_cl=n(5)./0.0355./n(1);
I=0.5*(m_na*(1)^2+m_cl*(-1)^2);
loggamma7=(ap1+ap2*Tk+ap3/Tk)*I-(ap4+ap5*Tk)*I/(I+1);
gamma7=exp(loggamma7);
elseif gammac=='duan_sun2003'
% Duan and Sun 2003 model; see ''Numerical modeling of geological carbon  
% sequestration: enhanced dissolution in randomly heterogeneous media (https://upcommons.upc.edu/handle/2117/376077)'' 
lamda7=-0.411370585+6.07632013e-4.*Tk+97.5347708./Tk-0.0237622469.*Pl./Tk+0.0170656236./(630-Tk).*Pl+1.41335834e-5.*Tk.*log(Pl);
xi7=3.36389723e-4-1.9829898e-5.*Tk+2.12220830e-3./Tk.*Pl-5.24873303e-3./(630-Tk).*Pl;
m_na=n(4)./0.023./n(1);
m_cl=n(5)./0.0355./n(1);
lngamma7=2.*lamda7.*m_na+xi7.*m_cl.*m_na;
gamma7=exp(lngamma7);
end
%% constant parameters

P0=1; %reference pressure 1 [bar]
R=83.144598; %universal gas constant []
bco2=27.8; % parameter to calculate the fugacity [cm^3/mol]
aco2=7.54e7-4.13e4*Tk; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
am=aco2; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
bm=bco2; % parameter to calculate the fugacity [cm^3/mol]
pv=[1, -R.*Tk./Pg, -(R.*Tk.*bm./Pg-am./Pg./Tk.^0.5+bm.^2), -am.*bm./Pg./Tk.^0.5];
v=roots(pv);
v=v(imag(v)==0);
v=sort(v);
vg=v(end);
vl=v(1);
w1=Pg.*(vg-vl);
w2=R.*Tk.*log((vg-bm)/(vl-bm))+am./Tk.^0.5./bm.*log((vg+bm).*vl./(vl+bm)./vg);
if w2<w1
    v=vl;
else
    v=vg;
end
logphi2=log(v/(v-bm))+bco2/(v-bm)+am*bco2/R/Tk^1.5/bm^2*(log((v+bco2)/v)...
-bm/(v+bm))-2*aco2/R/Tk^1.5/bm*log((v+bm)/v)-log(Pg*v/R/Tk);
% fugacity coefficient of gaseous CO2
phi2=exp(logphi2);
% fugacity of gaseous CO2
a2=phi2*Pg/P0;
% if disco2=0, dissolution is not considered and the dissolution is reduced
% by 1000 times


R_g=83.144598; %universal gas constant []
Tc=Tk-273.15;
if disco2==0
    df=1e3;
    Tc=40;
else
    df=1;
end
beta0=1.189+1.304e-2.*Tc-5.446e-5.*Tc.^2;
% equilibrium constant for CO2 evaporation;see  ''Numerical modeling of geological carbon  
% sequestration: enhanced dissolution in randomly heterogeneous media (https://upcommons.upc.edu/handle/2117/376077)'' 
k2=10^(beta0)*exp((Pg-1)*32.6/R_g/Tk)*df;
n(7)=n(1).*(a2./k2)./gamma7./55.5556;



end
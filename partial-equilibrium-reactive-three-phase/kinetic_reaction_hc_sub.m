function [N,u_m,r_k,dr_k_pg]=kinetic_reaction_hc_sub(Pl,Pg,Peq,nu,G,U,N,u,Tk,Sr,dt,rock_mass,Pg0,Ea1,Ea2,nr_chemical,disco2)


% SYNOPSIS:

% [N,u_m,r_k,dr_k_pg]=kinetic_reaction_hc_sub(Pl,Pg,Peq,nu,G,U,N,u,Tk,Sr,dt,rock_mass,Pg0,Ea1,Ea2,nr_chemical,disco2)
% DESCRIPTION:
%   This function serves to calculate the mass fraction of each species and components.
% the chemical system is
%    kinetic reaction  1:      CO2(g)+n_c*H2O(l)<-->HYD(s)
% equilibrium reaction 1:             CO2(aq)<==>CO2(g)
% equilibrium reaction 2:             NaCl(s)<==>Na+ + Cl-
% the stoichimetric matrix is
%         H2O(l) CO2(g) NaCl(s) Na+ Cl- HYD(s) CO2(aq)
%     K1:   -nc    -1     0      0   0    1      0
% nu= E1:   0      1     0      0   0    0      -1
%     E2:   0      0     -1     1   1    0      0
% The IDs of the   H2O(l), CO2(g), NaCl(s), Na+, Cl-, HYD(s) and CO2(aq) are
% 1 2 3 4 5 6 7, respectively,i.e.,  iH2O=1, iCO2g=2, iNaCl=3, iNa=4, iCl=5, iHYD=6 and iCO2aq=7 For instance, N(7) is the mass of CO2(aq).
% The kernel matrix is
%    H2O(l) CO2(g) NaCl(s) Na+ Cl- HYD(s) CO2(aq)
%  w  1      0      0       0   0   0.7105  0
%U=c  0      1      0       0   0   0.2895  0
%  Cl 0      0      1       0 1.6479   0    0
%  Z  0      0      0       1 -0.6479  0    0
% where w, c, Cl and Z are the four conservative components, of which the masses do not change during reaction.
% U satisfies U*MM*nu=0,
% where MM stores the molecular weight of all species, given as
% MM=diag([18,44,58.5,23,35.5,152,44].*1e-3).
% The governing equations are
% d n(iHYD)/d t= r_k
%     U*MM*n=u
% nu_e*ln(a)=ln(K_e)
%where, n(iHYD) is the molar abundance of hydrate, r_k is the hydration formation rate, u =[u_w,u_c,u_Cl,u_Z] stores the masses of the four components
%nu_e is the stoichiometric matrix for the two equilibrium reactions:
%           H2O(l) CO2(g) NaCl(s) Na+ Cl- HYD(s) CO2(aq)
% nu_e= E1:   0      1     0      0   0    0      -1
%       E2:   0      0     -1     1   1    0      0
% a stores the activities of all species, and K_e stores the two equilibrium
% constants for the two equilibrium reactions.
% NOMENCLATURE:
% Pl       -liquid pressure
% Pg       - gas pressure
% Peq      - equilibrium pressure for hydrate formation
% nu       -stoichiometric matrix
% G        - kernel matrix based on molar abundance
% U         -kernel matrix based on mass
% N         -mass of the species
% n         - molar abundance of the species
% u,u_m     - mass of the components
% Tk        -temperature in K
% Sr        - V*phi^1.5*(2.*kappa)^(-.5)*k0, with V being the volume, phi the porosity, kappa the intrinsic permeability and k0 the kinetic rate constant;
% dt        - time step
% rock_mass - the mass of the rock
% Pg0       - the gas pressure of previous time step
% Ea1, Ea2  - the activation energies for hydration
% nr_chemical -the threshold mass of hydration; when the water mass or the co2 mass is below nr_chemical is below nr_chemical, then hydration stops
% disco2      -indicator of whether we consider the CO2 dissolution; disco2=0, means we do not consider CO2 dissolution
% RETURNS:
% N   - mass of all species
% u_m - mass of all components
% r_k - hydrate formation rate
%dr_k_pg - derivative of hydrate formation rate with respect to gas
%pressure
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


%% stoichiometric matrices for kinetic and equilibrium reactions
nu_k=nu(1,:);
nu_e=nu(2:3,:);
%% change pressure from [pa] to [bar]
Pl=Pl.*1e-5;
Pg=Pg.*1e-5;
Peq=Peq.*1e-5;
Pg0=Pg0.*1e-5;
%% mass of species and component in [kg]
N=N';
u_m=u';
%% ID of seven species
iH2O=1; %iCO2g=2; iNaCl=3; iNa=4; iCl=5; iHYD=6; iCO2aq=7; 
% make sure the water mass is positive
if N(iH2O)<0
    N(iH2O)=1e-23;
end
%% diagnal matrix storing the molecular weight of all species 
MM=diag([18,44,58.5,23,35.5,152,44].*1e-3);
M=diag(1./([18,44,58.5,23,35.5,152,44].*1e-3));
%change from [kg] to [mol]
n=M*N;
n0=n;
% the molar volume of CO2 gas
v=VCO2S(Pg,Tk);

[A,a,loga]=aloga(n,v,Pg,Pl,Tk);
[r_k,dr_k,dr_k_pg]= rkdrk(nu_k,a,A,Tk,Pg0,Peq,n,Sr,dt,rock_mass,Ea1,Ea2,nr_chemical);

dtt=dt;
err=1;
nkk=1;
while err>1e-25&&nkk<20
    nkk=nkk+1;
    [L,R,salting]=LRH(dtt,n,n0,U,nu_e,nu_k,r_k,dr_k,A,MM,u_m,loga,Tk,Pg,disco2);   
    dn=-L\R;
    nold=n;
    if salting==0
        
        dn(3:5)=0;
        dn(5)=u(3)./0.0355/U(3,5)-nold(5);
        dn(4)=u(3)./0.0355/U(3,5)-nold(4);
        
    end
    
    ww=1e-1;
    eill=1e-8;
    theta=0.8;
    if norm(L*dn+R)/norm(-R)>eill
        dn=theta*dn;
        
        
    end
    
    indd=(dn+nold)<=0;
    dn(indd)=(ww-1).*nold(indd);
    n=dn+nold;
    err=max(abs(dn./(n+1e-0)));
    [A,a,loga]=aloga(n,v,Pg,Pl,Tk);
    [r_k,dr_k,dr_k_pg]= rkdrk(nu_k,a,A,Tk,Pg0,Peq,n,Sr,dt,rock_mass,Ea1,Ea2,nr_chemical);
    
end

N=MM*n;
u_m=U*N;
N=N';

end
function [A,a,loga]=aloga(n,v,Pg,Pl,Tk)
%% notation
%{
In the following, order of the species is H2O(l),CO2(g), NaCl(s), Na+, Cl-,
HYD(s) & CO2(aq). If we see  property(1), then it denotes the property of
H2O(l), and if we see property(2), then it denotes the property of
CO2(g), etc. For instance, m(7) denotes the molality of CO2(aq). 
%}

%% basic parameters
% number of species
ns=7;
%reference pressure 1 [bar]
P0=1; 
%universal gas constant []
R=83.144598; 
% parameter to calculate the CO2 fugacity [cm^3/mol] (c.f. Wang Yufei, 2022, Numerical modeling of geological carbon sequestration : enhanced dissolution in randomly heterogeneous media)
bco2=27.8; 
% parameter to calculate the CO2 fugacity [bar*cm^6*K^0.5/mol] (c.f. Wang Yufei, 2022, Numerical modeling of geological carbon sequestration : enhanced dissolution in randomly heterogeneous media)
aco2=7.54e7-4.13e4*Tk;
am=aco2; 
bm=bco2;
% molality of all species with respect to water
m=n./n(1).*55.5556;
 %molality of NaCl.BE CAREFUL nsm=0, but we should not use 0. and instead we should use 0.00001>0
nsm=0.5*(m(4)+m(5));
m_na=nsm;
m_cl=nsm;


%% calculate the activity of H2O(l) and the derivative of ln(activity): a1 and dlna1

a1=n(1)/sum(n([1 4 5 7]));
dlna1=zeros(7,1);
dlna1(1)=1/n(1)-1/sum(n([1 4 5 7]));
dlna1([ 4 5 7])=-1/sum(n([1 4 5 7]));
%% calculate the activity of CO2(g) and the derivative of ln(activity): a2 and dlna2
logphi2=log(v/(v-bm))+bco2/(v-bm)+am*bco2/R/Tk^1.5/bm^2*(log((v+bco2)/v)...
    -bm/(v+bm))-2*aco2/R/Tk^1.5/bm*log((v+bm)/v)-log(Pg*v/R/Tk);
% gas activity (or fugacity) coefficient
phi2=exp(logphi2);
% gas activity (or fugacity)
a2=phi2*Pg/P0;
dlna2=zeros(7,1);
%% calculate the activity of NaCl(s) and the derivative of ln(activity): a3 and dlna3  
% solid species is pure and activity is unit
a3=1;
dlna3=zeros(7,1);
%% calculate the activity of Na+ and the derivative of ln(activity): a4 and dlna4  
m_na(m_na==0)=1e-10;
m_cl(m_cl==0)=1e-10;
% we set activity coefficient =1; advanced model for activity can be found
% in Wang Yufei, 2022, Numerical modeling of geological carbon sequestration : enhanced dissolution in randomly heterogeneous media
a4=m(4);
dlna4=zeros(7,1);
dlna4(4)=dlna4(4)+1/n(4);
dlna4(1)=dlna4(1)-1/n(1);
%% calculate the activity of Cl- and the derivative of ln(activity): a5 and dlna5  
% we set activity coefficient =1
a5=m(5);
dlna5=zeros(7,1);
dlna5(5)=dlna5(5)+1/n(5);
dlna5(1)=dlna5(1)-1/n(1);
%% calculate the activity of hydrate and the derivative of ln(activity): a6 and dlna6 
% solid species is pure and activity is unit
a6=1;
dlna6=zeros(7,1);
%%  calculate the activity of CO2(aq) and the derivative of ln(activity): a7 and dlna7 
% the dissolution model is based on Duan and Sun 2003; another model based
% on Drummond (1981) has also been implelented in our previous versions. 
% c.f.  Wang Yufei, 2022, Numerical modeling of geological carbon sequestration : enhanced dissolution in randomly heterogeneous media
lamda7=-0.411370585+6.07632013e-4.*Tk+97.5347708./Tk-0.0237622469.*Pl./Tk+0.0170656236./(630-Tk).*Pl+1.41335834e-5.*Tk.*log(Pl);
xi7=3.36389723e-4-1.9829898e-5.*Tk+2.12220830e-3./Tk.*Pl-5.24873303e-3./(630-Tk).*Pl;
lngamma7=2.*lamda7.*m(4)+xi7.*m(5).*m(4);
gamma7=exp(lngamma7);
I=0.5*(m_na*(1)^2+m_cl*(-1)^2);
if I>5
    gamma7=5;
end
a7=gamma7.*m(7);
dlna7=zeros(ns,1);
dlna7(1)=-1./n(1); dlna7(7)=1./n(7);
%% 
A=[dlna1';dlna2';dlna3';dlna4';dlna5';dlna6';dlna7'];
a=[a1;a2;a3;a4;a5;a6;a7];
loga=log([a1;a2;a3;a4;a5;a6;a7]);
end
function [L,R,salting]= LRH(dt,n,n0,G_et,nu_e,nu_k,r_k,dr_k,A,MM,h,loga,Tk,Pg,disco2)


R_g=83.144598; %universal gas constant []
Tc=Tk-273.15;

if disco2==0
    df=1e3;
    Tc=40;
elseif  disco2==0.5 %geoxim
    df=1.08;
else
    df=1;
end

    
    
beta0=1.189+1.304e-2.*Tc-5.446e-5.*Tc.^2;
k2=10^(beta0)*exp((Pg-1)*32.6/R_g/Tk)*df;% equilibrium constant for co2 dissolution ;for Tc =40 We assume pl=pg
k3=10^(1.582);%equilibrium constant for salt dissolution
k=[k2;
    k3];% vector for equilibrium constant
L1=zeros(1,7);% generate space the first row of Jacobian matrix, i.e., the row for the hydrate formation
L1(6)=1/dt; %
L1=L1-(1).*dr_k';%
L=   [L1;
    G_et*MM;
    nu_e*A]; % Jacobian matrix
R=[(n(6)-n0(6))./dt-(1)*r_k;
    G_et*MM*n-h;
    nu_e*loga-log(k)]; % residual vector

a45=exp(loga([4,5]));
salting=1;
if a45(1)*a45(2)<k3 &&n(3)==0
    L(7,:)=[0 0 0 0 1 0 0];
    R(7)=0;
    salting=0;
end

end


function [r_k,dr_k,dr_k_pg]= rkdrk(nu_k,a,A,Tk,Pg,Peq,n,AS,dt,rock_mass,Ea1,Ea2,nr_chemical)
% the hydrate formation rate is calculated with
%r_k=-AS*Ar*Gamma_r*exp(Ea/(R*Tk))*(Peq-Pg)



%ID of seven speices
iH2O=1; iCO2g=2; iNaCl=3; iNa=4; iCl=5; iHYD=6; iCO2aq=7; 
R=8.314; %in j/K/mol; gas constant
if (n(iH2O)/6<nr_chemical||n(iCO2g)<nr_chemical) && Pg>Peq
    AS=0; % if the molar abundance is smaller than limit, stop reaction
end
if n(iHYD)<nr_chemical&&Pg<=Peq
    AS=0;% if molar abundance is smaller than limit, stop reaction
end
MH2O=0.018;MCO2=0.044;MHYD=0.152;% molecular weight of water, CO2 and hydrate
if Pg>Peq %hydrate formation, and the reaction is limited by water or CO2
    Ea=Ea1;
    if n(iH2O)/6<n(iCO2g)
        AS1=AS*(n(iH2O).*MH2O-nr_chemical.*6.*MH2O)./(rock_mass);% limited by water
    else
        AS1=AS*(n(iCO2g).*MCO2-nr_chemical.*MCO2)./(rock_mass);  %limited by gas
    end
else % hydrate dissociation and the reaction is limited  by HYD
   Ea=Ea2;
    AS1=AS.*(n(iHYD).*MHYD-nr_chemical.*MHYD)./(rock_mass);%limited by hydrate
    
end



r_k=AS1.*exp(-Ea./R./Tk).*(Pg-Peq).*1e5;%1e5 is for pressure form bar to pa
dr_k=zeros(7,1);
if Pg>Peq
    if n(iH2O)<n(iCO2g)
        dr_k(iH2O)=dr_k(iH2O)+MH2O./rock_mass*AS*exp(-Ea./R./Tk).*(Pg-Peq).*1e5;% limited by water
    else
        dr_k(iCO2g)=dr_k(iCO2g)+MCO2./rock_mass.*AS.*exp(-Ea./R./Tk).*(Pg-Peq).*1e5;%limited by gas
    end
else
    
    dr_k(iHYD)=dr_k(iHYD)+(MHYD)./(rock_mass).*...
        AS.*exp(-Ea./R./Tk).*(Pg-Peq).*1e5;% limited by hydrate
end

dr_k_pg=AS1*exp(-Ea./R./Tk)*1E5;
end

function v=VCO2S(Pg,Tk)
%% CO2 volume
% the gas volume is calculated based on the Redlich-Kwong equation, for detailed information, check the appendix B of  
%Spycher, N., Pruess, K., & Ennis-King, J. (2003). CO2-H2O mixtures in the geological sequestration of CO2. I. Assessment and calculation of mutual solubilities from 12 to 100 C and up to 600 bar. Geochimica et cosmochimica acta, 67(16), 3015-3031.
R=83.144598; %universal gas constant []
bco2=27.8; % parameter to calculate the fugacity [cm^3/mol]
aco2=7.54e7-4.13e4*Tk; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
am=aco2; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
bm=bco2; % parameter to calculate the fugacity [cm^3/mol]
% the governing  Redlich-Kwong equation is 
%v^3*1-v^2*(R*Tk/Pg)-v*(R*Tk*bm/Pg-am/Pg/Tk^0.5+bm^2) -am*bm/Pg/Tk^0.5=0
pv=[1, -R*Tk/Pg, -(R*Tk*bm/Pg-am/Pg/Tk^0.5+bm^2), -am*bm/Pg/Tk^0.5];
v=roots(pv);
v=v(imag(v)==0);
v=sort(v);
vg=v(end);
vl=v(1);
% if (w2 –w1)>0, then V is taken as the maximum root  and CO2 is in gas status; if (w2 – w1)<0, then
%V is taken as the minimum root and CO2 is in liquid status.
w1=Pg*(vg-vl);
w2=R*Tk*log((vg-bm)/(vl-bm))+am/Tk^0.5/bm*log((vg+bm)*vl/(vl+bm)/vg);
if w2<w1
    v=vl;
else
    v=vg;
end
end
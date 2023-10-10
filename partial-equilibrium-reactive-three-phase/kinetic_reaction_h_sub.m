 function [n,h,r_k,dr_k_pg]=kinetic_reaction_h_sub(Pl,Pg,Peq,nu,G,G_t,n,H,Tk,Sr,dt,rock_mass,ksal,Pg0,disco2)
 
 % speciation calculation based on the partial equilibrium reaction given in section called
%Partial Equilibrium Reactive Multi-Component Three-Phase Flow Mode in the user'guides
%
% SYNOPSIS:
%[n,h]=kinetic_reaction_sub(Pl,Pg,nu,G,G_t,n,H,Tk,Sr,dt)
%
% DESCRIPTION:
%   This function serves to calculate the mass fraction of each species and components. 
%Description is already included in the user's guide. 

%
% REQUIRED PARAMETERS:
% Pl       -liquid pressure
%Pg      - gas pressure
% nu   -stoichiometric matrix
% G        - kernel matrix based on molar abundance
% G_t      -kernel matrix based on mass
% n         -mass of the species
% H         - mass of the components
% Tk       -temperature in K
%Sr, dt   - useless for equilibrium reactions

% RETURNS:
%   n, h   -  mass of species and components

% SEE ALSO:
%  kinetic_reaction

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
 
 
%% born coefficient and effective electrostatic radius

nu_k=nu(1,:);
nu_e=nu(2:3,:);
%%


Pl=Pl.*1e-5;
Pg=Pg.*1e-5;
Peq=Peq.*1e-5;
Pg0=Pg0.*1e-5;
%% Right hand side


n=n';
h=H';

%n(1)=(h(1)-G_t(1,6)*n(6))/G_t(1,1);
if n(1)<0
    n(1)=1e-23;
end


nk=1;
MM=diag([18,44,58.5,23,35.5,152,44].*1e-3);
M=diag(1./([18,44,58.5,23,35.5,152,44].*1e-3));
n=M*n;
n0=n;
t=0;
v=VCO2S(Pg,Tk);
 [A,a,loga]=aloga(n,v,Pg,Tk);
 [r_k,dr_k,dr_k_pg]= rkdrk(nu_k,a,A,Tk,Pg0,Peq,n,Sr,dt,rock_mass);
 


phase_3=2;
if phase_3==1% equilibrium three phase
n=reaction_kinetic_equilibrium(Pg,Pg,nu,G_t,n,H,Tk); 
end
   warning off

%%
if phase_3==2% true kinetic

 dtt=dt;

 while t<dt
     err=1;
  nkk=1;
while err>1e-25&&nkk<12
nkk=nkk+1;
[L,R,salting]=LRH(dtt,n,n0,G_t,nu_e,nu_k,r_k,dr_k,A,MM,h,loga,Tk,Pg,disco2);


dn=-L\R;

if salting==0
    
   dn(3:5)=0;   

end
nold=n;


ww=eps;
ww=0.1;
eill=1e-8;
theta=1;
if norm(L*dn+R)/norm(-R)>eill
    dn=theta*dn;
    
    
end

indd=(dn+nold)<=0;
dn(indd)=(ww-1).*nold(indd);
n=dn+nold;   
err=max(abs(dn./(n+1e-20)));
    [A,a,loga]=aloga(n,v,Pg,Tk);
   [r_k,dr_k,dr_k_pg]= rkdrk(nu_k,a,A,Tk,Pg0,Peq,n,Sr,dt,rock_mass);

end

nk=nk+1;
t=t+dtt;


 dtt=max(dtt,dt/5);
 dtt=min(dtt,dt-t);
  % r_k*dtt*6>n(1)
 
 end
 %%
%  if isnan(n(1))||nkk==13
%      t=0; dtt=dt/10;
%      n=n0;
%       [A,a,loga]=aloga(n,v,Pg,Tk);
%    [r_k,dr_k]= rkdrk(nu_k,a,A,Tk,Pg,Peq,n,Sr,dt,rock_mass);
%   while t<dt
%      err=1;
%   nkk=1;
% while err>1e-12&&nkk<12
% nkk=nkk+1;
% [L,R,salting]=LRH(dtt,n,n0,G_t,nu_e,nu_k,r_k,dr_k,A,MM,h,loga,Tk,Pl);
% 
% 
% dn=-L\R;
% 
% if salting==0
%     
%     dn(3:5)=0;   
% 
% end
% nold=n;
% 
% 
% ww=0.6;
% eill=1e-8;
% theta=1;
% if norm(L*dn+R)/norm(-R)>eill
%     dn=theta*dn;
%     
%     
% end
% 
% indd=(dn+nold)<=0;
% dn(indd)=(ww-1).*nold(indd);
% n=dn+nold;   
% err=max(abs(dn./(n+1e-20)));
%     [A,a,loga]=aloga(n,v,Pg,Tk);
%    [r_k,dr_k]= rkdrk(nu_k,a,A,Tk,Pg,Peq,n,Sr,dt,rock_mass);
% 
% end
% 
% nk=nk+1;
% t=t+dtt;
% 
% 
%  dtt=min(dtt,dt-t);
%   
%  n0=n;
%  end
%  
%  end
 

   % ntt=n;
  % n=reaction_kinetic_equilibrium(Pl,Pg,nu,G_t,n,H,Tk); 
    

end
n=MM*n;
h=G_t*n;
n=n';

end
function [A,a,loga]=aloga(n,v,Pg,Tk)
%% born coefficient and effective electrostatic radius
m=n./n(1).*55.5556;
r_na=1.91;r_cl=1.81;r_h=3.08;r_ca=2.87; r_hco3=2.10;r_co3=2.81;r_oh=1.4;
w_na=0.3305*10^5;w_cl=1.456*10^5;w_h=0.5387*10^5;w_ca=2.314*10^5;w_hco3=0.7906*10^5;w_co3=2.3634*10^5;w_oh=1.1859*10^5;
ns=7;
z5=2;z6=-1;z8=-2;z9=-1;
%%

P0=1; %reference pressure 1 [bar]

Tc=Tk-273.15;%temperature in [kelvin]
R=83.144598; %universal gas constant []
bh2o=18.18; % parameter to calculate the fugacity [cm^3/mol]
bco2=27.8; % parameter to calculate the fugacity [cm^3/mol]
ach=7.89e7; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
aco2=7.54e7-4.13e4*Tk; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
am=aco2; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
bm=bco2; % parameter to calculate the fugacity [cm^3/mol]
%% calculate a1 and dlnaa;H2O(l)
nsm=0.5*(m(4)+m(5)); %molality of NaCl.BE CAREFUL nsm=0, but we should not use 0. and instead we should use 0.00001>0
nsn=[n(1)*18.0/1000*nsm];%
m_na=nsm;
m_cl=nsm;
n_na=nsn;
n_cl=nsn;
I_m=0.5*(m_na*(1)^2+m_cl*(-1)^2);
I_m(I_m==0)=1e-10;
dI_m=zeros(7,1);
dI_m(4)=0.5*1*55.5556./n(1);dI_m(5)=0.5*1*55.5556/n(1);
dI_m(1)=-(dI_m(4)*n(4)+dI_m(5)*n(5))/n(1);
x_w=n(1)/(n(1)+n(4)+n(5)+n(7));% should we consider the mole of ion. it does not matter because the number of ion is small
dx_w=zeros(ns,1);
dx_w(1)=(n(4)+n(5)+n(7))./(n(4)+n(5)+n(7)+n(1)).^2;
dx_w(4)=-n(1)./(n(4)+n(5)+n(7)+n(1)).^2;dx_w(5)=-n(1)./(n(4)+n(5)+n(7)+n(1)).^2;dx_w(7)=-n(1)./(n(4)+n(5)+n(7)+n(1)).^2;
%a_gamma=input('from table B1 of allan');
a_gamma=0.5281;
%b_nacl=input('from table B.3 of allan');
b_nacl=4.783e-7;
%b_na_cl=input('from table B.4 of allan');
b_na_cl=-5.603e-2;
%b_gamma=input('from table B2 of allan');
b_gamma=0.5281e8;



a_o=(m_na*r_na+m_cl*r_cl)/(m_na+m_cl);
AA=1+a_o*b_gamma*(I_m)^0.5;
sigma=3/(a_o*b_gamma*(I_m)^0.5)^3*(AA-1/AA-2*log(AA));
psi_na=a_gamma*(1)^2*(I_m)^0.5/3*sigma+x_w/(1-x_w)*log10(x_w)-0.5*(w_na*b_nacl+b_na_cl-0.19*(1-1))*I_m;
psi_cl=a_gamma*(-1)^2*(I_m)^0.5/3*sigma+x_w/(1-x_w)*log10(x_w)-0.5*(w_cl*b_nacl+b_na_cl-0.19*(1-1))*I_m;
loga1=2.303/55.508*(m_na*psi_na+m_cl*psi_cl);
a1=exp(loga1);
%dx_w_2=(n(4)+n_na+n_cl)/(n(2)+n(4)+n_na+n_cl)^2;
%dx_w_4=-n(2)/(n(2)+n(4)+n_na+n_cl)^2;
dx_w_2=(n(7)+n_na+n_cl)/(n(1)+n(7)+n_na+n_cl)^2;
dx_w_4=-n(1)/(n(1)+n(7)+n_na+n_cl)^2;
dlna1_2=2.303/55.508*((dx_w_2/(1-x_w)^2*log10(x_w)+x_w/(1-x_w)*log10(2.7183)/x_w*dx_w_2)*m_na+...
    (dx_w_2/(1-x_w)^2*log10(x_w)+x_w/(1-x_w)*log10(2.7183)/x_w*dx_w_2)*m_cl);

dlna1_4=2.303/55.508*((dx_w_4/(1-x_w)^2*log10(x_w)+x_w/(1-x_w)*log10(2.7183)/x_w*dx_w_4)*m_na+...
    (dx_w_4/(1-x_w)^2*log10(x_w)+x_w/(1-x_w)*log10(2.7183)/x_w*dx_w_4)*m_cl);

dlna1=zeros(ns,1);
dlna1(1)=dlna1_2; dlna1(7)=dlna1_4; 

a1=n(1)/sum(n([1 4 5 7]));dlna1=zeros(7,1);
dlna1(1)=1/n(1)-1/sum(n([1 4 5 7]));
dlna1([ 4 5 7])=-1/sum(n([1 4 5 7]));
 %% calculate a2 and dlna2;CO2(g)
 logphi2=log(v/(v-bm))+bco2/(v-bm)+am*bco2/R/Tk^1.5/bm^2*(log((v+bco2)/v)...
-bm/(v+bm))-2*aco2/R/Tk^1.5/bm*log((v+bm)/v)-log(Pg*v/R/Tk);

phi2=exp(logphi2);

a2=phi2*Pg/P0;
dlna2=zeros(7,1);

%% calculate a3 and dloga3;NaCl(s)
a3=1;
dlna3=zeros(7,1);
%% calculate a4 and dlna4;Na+
m_na(m_na==0)=1e-10;
m_cl(m_cl==0)=1e-10;
a_o=(m_na*r_na+m_cl*r_cl)/(m_na+m_cl);
AA=1+a_o*b_gamma*(I_m)^0.5;

% log10gamma4=-a_gamma*(1)^2*sqrt(I_m)/AA+log10(x_w)+(w_na*b_nacl+b_na_cl-0.19*(abs(1)-1))*I_m;
% a4=m(4)*10^(log10gamma4);
% d_AA=a_o.*b_gamma.*0.5./sqrt(I_m).*dI_m;
% dlna4=log(10)*((-a_gamma*(1)^2).*(0.5./sqrt(I_m).*dI_m./AA-d_AA.*sqrt(I_m))./AA.^2+1./log(10)./x_w.*dx_w+(w_na*b_nacl+b_na_cl-0.19*(abs(1)-1))*dI_m);
%  if 10^(log10gamma4)>3
%  a4=m(4)*3.0;
% dlna4=zeros(7,1);
%  end
% 
% dlna4(4)=dlna4(4)+1/n(4);
% dlna4(1)=dlna4(1)-1/n(1);
%% if we set activity coefficient =1
a4=m(4);
dlna4=zeros(7,1);
dlna4(4)=dlna4(4)+1/n(4);
dlna4(1)=dlna4(1)-1/n(1);
%% calculate a5 and dlna5;cl-
% log10gamma5=-a_gamma*(-1)^2*sqrt(I_m)/AA+log10(x_w)+(w_cl*b_nacl+b_na_cl-0.19*(abs(-1)-1))*I_m;
% a5=m(5)*10^(log10gamma5);
% dlna5=log(10)*((-a_gamma*(-1)^2).*(0.5./sqrt(I_m).*dI_m./AA-d_AA.*sqrt(I_m))./AA.^2+1./log(10)./x_w.*dx_w+(w_cl*b_nacl+b_na_cl-0.19*(abs(-1)-1))*dI_m);
%  if 10^(log10gamma5)>3
%  a5=m(5)*3.0;
% dlna5=zeros(7,1);
%  end
% dlna5(5)=dlna5(5)+1/n(5);
% dlna5(1)=dlna5(1)-1/n(1);
%%
a5=m(5);
dlna5=zeros(7,1);
dlna5(5)=dlna5(5)+1/n(5);
dlna5(1)=dlna5(1)-1/n(1);
%% calculate a6 and dlna6;HYD
a6=1;
dlna6=zeros(7,1);
%% calculate a7 and dlna7;co2(aq)
ap1=-1.0312;
ap2=1.2806e-3;
ap3=255.9;
ap4=0.4445;
ap5=-1.606e-3;
I=0.5*(m_na*(1)^2+m_cl*(-1)^2);
dI1=-0.5.*(n(4)./n(1)^2./.018+n(5)./n(1)^2./.018);
dI45=0.5./n(1)./0.018;
loggamma7=(ap1+ap2*Tk+ap3/Tk)*I-(ap4+ap5*Tk)*I/(I+1);
gamma7=exp(loggamma7);
if I>5
    gamma7=5;dI1=0;dI45=0;
end
a7=55.5556*n(7)/n(1)*gamma7;
dlna7=zeros(ns,1);
dlna7(1)=-1./n(1); dlna7(7)=1./n(7);

dlngamma7_1=((ap1+ap2*Tk+ap3/Tk).*dI1-(ap4+ap5*Tk)./(I+1).^2.*dI1);
dlngamma7_45=((ap1+ap2*Tk+ap3/Tk).*dI45-(ap4+ap5*Tk)./(I+1).^2.*dI45);
dlna7(1)=dlna7(1)+dlngamma7_1;
dlna7(4)=dlngamma7_45;
dlna7(5)=dlngamma7_45;
%%

A=[dlna1';dlna2';dlna3';dlna4';dlna5';dlna6';dlna7'];
a=[a1;a2;a3;a4;a5;a6;a7];
loga=log([a1;a2;a3;a4;a5;a6;a7]);
% id=(n<10^(-16))&(n>0);
% A(id,id)=1./n(id);


end
function [L,R,salting]= LRH(dt,n,n0,G_et,nu_e,nu_k,r_k,dr_k,A,MM,h,loga,Tk,Pg,disco2)


R_g=83.144598; %universal gas constant []
if disco2==0
k2=10^(1.623)*exp((Pg-1)*32.6/R_g/Tk)*1e3;%for Tc =40 We assume pl=pg
else
k2=10^(1.623)*exp((Pg-1)*32.6/R_g/Tk);
end

k3=10^(1.582);

% if n(3)==0
%     k3=a45(1)*a45(2);
% end
k=[k2;
    k3];
L1=zeros(1,7);
% L1(6)=1/dt;
% L1=L1-dr_k';

L1(6)=1/dt;
L1=L1-(1).*dr_k';

L=   [L1;
      G_et*MM;
     nu_e*A];
 
 

% R=[(n(6)-n0(6))./dt-(1)*r_k;
%      G_et*MM*n-h;
%      nu_e*loga-log(k)];
 
 R=[(n(6)-n0(6))./dt-(1)*r_k;
     G_et*MM*n-h;
     nu_e*loga-log(k)];
 
 a45=exp(loga([4,5]));
salting=1;
if a45(1)*a45(2)<k3 &&n(3)==0
L(7,:)=[0 0 0 0 1 0 0];
R(7)=0;
   salting=0;
end

end


function [r_k,dr_k,dr_k_pg]= rkdrk(nu_k,a,A,Tk,Pg,Peq,n,AS,dt,rock_mass)


[Ea1,Ea2]=deal(0,0);
R=8.314; %in j/K/mol

if (n(1)*0.018/(n(1)*0.018+n(2)*0.044)<1e-5||n(2)*0.044/(n(1)*0.018+n(2)*0.044)<1e-5) && Pg>Peq
   AS=0; 
end

if Pg>Peq
%     AS1=AS.*n(1).*n(2)./(n(1)+n(2)).^2;
    
    if n(1)<n(2)
        AS1=AS*n(1).*0.018./(rock_mass);
    else
       AS1=AS*n(2).*0.044./(rock_mass);  
    end
    
    
else
    AS1=AS.*n(6).*0.152./(n(6).*0.152+n(3)*0.0585+rock_mass);
    
end

% AS1(n(6)<=0&(Pg-Peq)<0)=0;
% AS1((n(1)<=0|n(2)<=0)&(Pg-Peq)>0)=0;


r_k=AS1.*exp(-Ea1./R./Tk).*(Pg-Peq).*1e5;%1e5 is for pressure form bar to pa
% if r_k==0
%     if n(6)>0&&(Pg-Peq)<0
%     r_k=-n(6)/dt*1.5;
%     end
%     
%     if (n(1)>0&&n(2)>0)&&(Pg-Peq)>0
%         r_k=min(n(1)/6,n(2))/dt*1.5;
%     end
% end

dr_k=zeros(7,1);
if Pg>Peq
    
  
%     dr_k(1)=dr_k(1)+(n(2).^3-n(1).^2.*n(2))./(n(2)+n(1)).^4.*AS*exp(-Ea1./R./Tk).*(Pg-Peq).*1e5;
%       dr_k(2)=dr_k(2)+(n(1).^3-n(2).^2.*n(1))./(n(2)+n(1)).^4.*AS.*exp(-Ea1./R./Tk).*(Pg-Peq).*1e5;
      
      if n(1)<n(2)
           dr_k(1)=dr_k(1)+0.018./rock_mass*AS*exp(-Ea1./R./Tk).*(Pg-Peq).*1e5;
         
           
      else
      dr_k(2)=dr_k(2)+0.044./rock_mass.*AS.*exp(-Ea1./R./Tk).*(Pg-Peq).*1e5;
          
      end
      
    
else
   dr_k(3)=dr_k(3)+(-n(6)*0.0585*0.152)./(n(6).*0.152+n(3)*0.0585+rock_mass).^2.*...
        AS.*exp(-Ea1./R./Tk).*(Pg-Peq).*1e5;
    
    dr_k(6)=dr_k(6)+((n(6).*0.152+n(3)*0.0585+rock_mass).*0.152-0.152*0.152*n(6))./(n(6).*0.152+n(3)*0.0585+rock_mass).^2.*...
        AS.*exp(-Ea1./R./Tk).*(Pg-Peq).*1e5;
    
    
end
% if n(6)<=rock_mass*1e-4&&(Pg-Peq)<0
% r_k=0;
% 
% dr_k=zeros(7,1);
% 
% end
% if (n(1)<=rock_mass*1e-5||n(2)<=rock_mass*1e-5)&&(Pg-Peq)>0;
% r_k=0;
% dr_k=zeros(7,1);
% end

  dr_k_pg=AS1*exp(-Ea1./R./Tk)*1E5;
end

function v=VCO2S(Pg,Tk)
%% CO2 volume
P0=1; %reference pressure 1 [bar]
R=83.144598; %universal gas constant []
bh2o=18.18; % parameter to calculate the fugacity [cm^3/mol]
bco2=27.8; % parameter to calculate the fugacity [cm^3/mol]
ach=7.89e7; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
aco2=7.54e7-4.13e4*Tk; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
am=aco2; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
bm=bco2; % parameter to calculate the fugacity [cm^3/mol]
    pv=[1, -R*Tk/Pg, -(R*Tk*bm/Pg-am/Pg/Tk^0.5+bm^2), -am*bm/Pg/Tk^0.5];
v=roots(pv);
v=v(imag(v)==0);
v=sort(v);
vg=v(end);
vl=v(1);
w1=Pg*(vg-vl);
w2=R*Tk*log((vg-bm)/(vl-bm))+am/Tk^0.5/bm*log((vg+bm)*vl/(vl+bm)/vg);
if w2<w1
    v=vl;
else
    v=vg;
end
end

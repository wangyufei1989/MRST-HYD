function [d_rho,d_mu]=D_RHOMU_kinetic_h(state)
% Calculate the  density gradient and viscosity gradient of brine and co2 phase at
% given pressure, mass fraction, NaCl concentration and temperature.
%
% SYNOPSIS:
%  [d_rho,d_mu]=D_RHOMU_kinetic(state)
%
% DESCRIPTION:
% This function serves to calculate the  density gradient and viscosity gradient of brine and co2 phase at
% given pressure, mass fraction, NaCl concentration and temperature.
%The following model for brine denisty can be found in Garcia [2001;2003].
%the following model for brine viscosity can be found in Garcia [2003].
%the following model for CO2 density can be found in Victor Vilarrasa [2012] phd thesis
%or Spycher et al. [2003].
%the following model for CO2 viscosity can be found in Victor Vilarrasa [2012] phd thesis.
%Description is already included in the user's guide.
%
% REQUIRED PARAMETERS:
% state   -state varibles
% RETURNS:
%   d_rho   -density gradient of brine and gas over pressure- [drho_l/dp_l,drho_g/dp_g,drho_l/dp_g,drho_g/dp_l];
%
%   d_mu    -  viscosity gradient of brine and gas over pressure-[dmu_l/dp_l,dmu_g/dp_g,dmu_l/dp_g,dmu_g/dp_l];  .
%
% SEE ALSO:
%   rhomu_p_frac_kinetic_h

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
% the brine density derivative

if size(state.pressure,2)==1
    state.pressure=[state.pressure,state.pressure];
end
p0=state.pressure(:,1);
m_NaCl=state.m_NaCl;
m_NaCl=max(m_NaCl,1e-10);
mc=state.species(:,7)./(sum(state.species,2)-state.species(:,2)-state.species(:,3)-state.species(:,6));%mass fraction of co2
dsp=(state.dsp);
dsp=zeros(size(dsp));
dcp=(state.dcp);
% dsp=0;
% dcp=0;
Mc=0.044; %molar mass of co2
T_c=state.Tk-273.15;
Tk=T_c+273.15;
vphi=3.751*10^(-5)-9.585*10^(-8).* T_c+8.74*10^(-10)* T_c.^2-5.044*10^(-13).* T_c.^3;
x=-9.9595.*exp( -4.539*10^(-3).*m_NaCl )+7.0845*exp( -1.638*10^(-4).*T_c)...
    + 3.9093*exp(2.551*10^(-10).*p0);
rhor=-3.033405*10^3+1.0128163*10^4 *x -8.750567*10^3* x.^2+2.66310*10^3* x.^3;
dx_p=3.9093*exp(2.551*10^(-10).*p0).*2.551*10^(-10);
drhor_p=1.0128163*10^4 *dx_p-2*8.750567*10^3.* x.*dx_p+3*2.66310*10^3* x.^2.*dx_p;
% derivitive of water density with respect to the liquid pressure 
drhob_p=drhor_p.*(1+mc.*(1-rhor.*vphi./Mc))+rhor.*(-mc.*drhor_p.*vphi./Mc);
drhob_pg=rhor.*(1-rhor.*vphi./Mc).*dcp;
dx_salinity=-9.9595 *exp( -4.539*10^(-3)*m_NaCl ).*(-4.539*10^(-3)).*dsp;
drhor_salinity=1.0128163*10^4 *dx_salinity-2*8.750567*10^3.* x.*dx_salinity+3*2.66310*10^3* x.^2.*dx_salinity;
% derivitive of water density with respect to the gas pressure 
drhob_pg=drhob_pg+drhor_salinity.*(1+mc.*(1-rhor.*vphi./Mc))-rhor.*mc.*vphi./Mc.*drhor_salinity;
% Temperature
dx_t=7.0845.*exp(-1.638e-4.*T_c).*(-1.638e-4);
drhor_t=1.0128163*10^4 *dx_t-2*8.750567*10^3.* x.*dx_t+3*2.66310*10^3* x.^2.*dx_t;
dvphi_t=-9.585*10^(-8)+8.74*10^(-10).* T_c.*2-5.044*10^(-13).* T_c.^2.*3;
% derivitive of water density with respect to the temperature 
drhob_t=drhor_t.*(1+mc.*(1-rhor.*vphi./Mc))+rhor.*mc.*(-drhor_t.*vphi./Mc-rhor.*dvphi_t./Mc);
% now it is the co2 density derivative
p0=state.pressure(:,2);
[v,dv]=VCO2(p0,Tk);
% derivitive of CO2 density with respect to the gas pressure 
drhoc_p=-(0.044)./v.^2.*dv;
% derivitive of CO2 density with respect to the liquid pressure 
drhoc_pl=zeros(size(drhoc_p));
Tkd=(1-1e-7).*Tk;
[vd,~]=VCO2(p0,Tkd);
drhoc_t=(0.044./vd-0.044./v)./(-1e-7.*Tk);
%  [v,dvt]=VCO2t(p0,Tk);
%  drhoc_t=-(0.044)./v.^2.*dvt;

%derivative of viscosity of water
Mc=0.044; %molar mass of co2
m_C=state.species(:,7)./Mc./state.species(:,1); % molarity of co2
% derivitive of water viscosity with respect to the liquid pressure 
dMub_p=1.69956*10^(-9)*10^(-3).*ones(size(state.pressure,1),1);
% derivitive of water viscosity with respect to the gas pressure 
dMub_pg=((8.79552-3.17229*10^(-2)* Tk)+(-7.22796+2.64498*10^(-2) *Tk).* m_C.*2).*1e-3./0.044./(1-state.salinity).*dcp...
    +((3.85971-1.32561*10^(-2).* Tk).*dsp+(-5.37539+1.90621*10^(-2).* Tk).*...
    0.5.*m_NaCl.^(-0.5).*dsp).*1e-3;%approximation is used
% derivitive of water viscosity with respect to the temperature 
dMub_t=(-1.32561*10^(-2).* m_NaCl+(1.90621*10^(-2)).*...
    m_NaCl.^(1/2)+(-3.17229*10^(-2)).* m_C+(2.64498*10^(-2)).* m_C.^2).*10^(-3);

%viscosity of co2
tc=Tk/304;%temperature scaled by critical temperature.
a10  = 0.248566120;
a11  = 0.004894942;
a20  = -0.373300660;
a21  = 1.22753488;
a30  = 0.363854523;
a31 = -0.774229021;
a40  = -0.0639070755;
a41  = 0.142507049;
rhoc=state.rho(:,2)./468;%scaled by critical density
drhoc=drhoc_p./468;
% derivitive of CO2 viscosity with respect to the gas pressure
dMuc_p=state.mu(:,2).*(a10.*drhoc+a11*drhoc./tc+2.*a20.*rhoc.*drhoc+2.*a21*rhoc.*drhoc./tc...
    +3.*a30.*rhoc.^2.*drhoc+3.*a31*rhoc.^2.*drhoc./tc+4.*a40*rhoc.^3.*drhoc+4*a41*rhoc.^3.*drhoc./tc);
% derivitive of CO2 viscosity with respect to the liquid pressure
dMuc_pl=zeros(size(dMuc_p));
% derivitive of CO2 density with respect to the temperature 
drhoct=drhoc_t./468;
% derivitive of CO2 viscosity with respect to the temperature 
dMuc_t=state.mu(:,2).*(a10.*drhoct+a11.*drhoct./tc+2.*a20.*rhoc.*drhoct+2.*a21*rhoc.*drhoct./tc...
    +3.*a30.*rhoc.^2.*drhoct+3.*a31.*rhoc.^2.*drhoct./tc+4.*a40.*rhoc.^3.*drhoct+4.*a41.*rhoc.^3.*drhoct./tc...
    -a11.*rhoc./tc.^2./304-a21.*rhoc.^2./tc.^2./304-a31.*rhoc.^3./tc.^2./304-a41*rhoc.^4./tc.^2./304)...
    ...
    +((tc).^0.5.*(16.6346068./(tc).^2./304 -2.*4.66920556./(tc ).^3./304).*10^(-6)...
    +0.5.*(tc).^(-0.5)./304.*(27.2246461-16.6346068./(tc) +4.66920556./(tc ).^2).*10^(-6)).*exp(a10*rhoc+a11*rhoc./tc+a20*rhoc.^2+a21*rhoc.^2./tc...
    +a30*rhoc.^3+a31*rhoc.^3./tc+a40*rhoc.^4+a41*rhoc.^4./tc);

d_mu=[dMub_p,dMuc_p,dMub_pg,dMuc_pl,dMub_t,dMuc_t];
d_rho=[drhob_p,drhoc_p,drhob_pg,drhoc_pl,drhob_t,drhoc_t];
if isfield(state,'rhomu')
    if state.rhomu.rhomu_flag==1% use simplified density viscosity model
        drhoc_p=state.rhomu.rhor(2).*state.rhomu.compress(2).*ones(size(p0)).*1;
        drhoc_t=-state.rhomu.rhor(2).*state.rhomu.compress(4).*ones(size(p0)).*1;
        drhob_t=-state.rhomu.rhor(1).*state.rhomu.compress(3).*ones(size(p0)).*1;
        drhoc_pl=zeros(size(drhoc_p));
        drhob_pg=zeros(size(drhoc_p));
        drhob_p=state.rhomu.rhor(1).*state.rhomu.compress(1).*ones(size(p0)).*1;
        d_rho=[drhob_p,drhoc_p,drhob_pg,drhoc_pl,drhob_t,drhoc_t];
        d_mu=zeros(size(drhoc_p,1),6);
    end
    
end



end
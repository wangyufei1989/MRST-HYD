function fluid = initFluid_VG(varargin)


% SYNOPSIS:
% fluid = initFluid_VG(varargin)
%
% DESCRIPTION:

% Define the relative permeability and retention curves based on REVISED Van Genuchten model;
% define the hydration equilibrium pressure based on temperature and NaCl concentration.

%---------- relative permeability model------------%
% In this model the relative permeability is based on the simplest Corey model: 
% krl=k_lmax*((s_l-s_lr)/(s_lmax-s_lr))^n_l
% krg=k_gmax*((s_g-s_gr)/(s_gmax-s_gr))^n_g
% where krl= water relative permeability, k_lmax= maximum water relative
% permeability, s_l= water saturation, s_lr=residual water saturation,
% s_lmax= maximum water saturation, n_l= shape factor=2;
% similalry, krg= CO2 relative permeability, k_gmax= maximum CO2 relative
% permeability, s_g= CO2 saturation, s_gr=residual CO2 saturation,
% s_gmax= maximum CO2 saturation, n_g= shape factor=2;

%-----------------retention curve----------------%
% In this model the retention curve is based on the Van Genuchtten model: 

%     / 1, pc<0
%s_l=
%     \ 1/(1+(sqrt(kappa*phi_m/(kappa_g*phi))*alpha_van*pc)^np)^mp, pc>=0
%where pc=p_g-p_l=gas pressure - liquid pressure; kappa is the intrinsic
%permeability; phi is the porosity; kappa_g is the geometric mean
%permeability; phi_m is mean porosity; alpha_van is the scaling pressure
%np=1/(1-mp) is the shape factor

% REQUIRED PARAMETERS:
%'sr=[s_lr,s_gr]', residual water and CO2 saturations,
%'kwm=[k_lmax,k_gmax]', maximum water and CO2 relative permeability,
%'pc_scale', scaling parameter, 'alpham', shape factor ;


% RETURNS:
%   fluid   - an object containing relative permeability, retention curves and equilibrium hydration pressure.
%




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


opt = struct( 'sr', [], 'kwm',[], 'pc_scale',[],'alpham',[],'pcmax',[]);
opt = merge_options(opt, varargin{:});
kr  = @(s,varargin) relperm(s, opt, varargin{:});
S   = @(state) S_funct(state,opt);
peq=@(state) peq_funct(state);
fluid = struct('relperm'   , kr,...
                'S'        , S,...
                'peq'      ,peq);
end


%--------------------------retention curve-----------------------------------------------


function varargout = S_funct(state,opt)
scale=2;
opt.pcmax=opt.pcmax/(1+scale);
pc=(state.pressure(:,2)-state.pressure(:,1));
pcreal=pc;
pc(pc<0)=0;
pc_cut=pc>opt.pcmax;
pc(pc_cut)=opt.pcmax;
s_cut=zeros(size(pc));
p_d=0;


ps=opt.pc_scale;
ps=ps./10^5.*opt.alpham(1);
m=opt.alpham(2); n=1/(1-m);


opt.sr(1)=0;
s_cut(pc_cut)=(1+(ps(pc_cut).*opt.pcmax).^n).^(-m).*(1-opt.sr(1))+opt.sr(1)-ps(pc_cut).*p_d.*opt.pcmax;

varargout{1}   =[(1+(ps.*pc).^n).^(-m).*(1-opt.sr(1))+opt.sr(1)-ps.*p_d.*pc,...
    1-((1+(ps.*pc).^n).^(-m).*(1-opt.sr(1))+opt.sr(1)-ps.*p_d.*pc)];

pc_cut2=pcreal>opt.pcmax*(1+scale);
varargout{1}(pc_cut,1)= varargout{1}(pc_cut)-(pcreal(pc_cut)-opt.pcmax)./opt.pcmax.*s_cut(pc_cut)./scale;
varargout{1}(pc_cut,2)=1-( varargout{1}(pc_cut));
varargout{1}(pc_cut2,1)=0;
varargout{1}(pc_cut2,2)=1;
if nargout > 1, varargout{2} = -m.*(1+(ps.*pc).^n).^(-m-1).*ps.^n.*n.*pc.^(n-1).*(1-opt.sr(1))-ps.*p_d;
    
    varargout{2}(pc_cut)=-s_cut(pc_cut)./opt.pcmax./scale;
    varargout{2}(pc_cut2)=0;
end

end


%-----------------------equilibrium hydration pressure--------------%


function varargout = peq_funct(state)

Tk=state.Tk;
salinity=state.m_NaCl;
t=[200 268.15 269.15 270.15 271.15 272.15 273.15 274.15 275.15 276.15 277.15 278.15 279.15 280.15 281.15 282.15 283.15 284.15 285.15 286.15 287.15 288.15 310 350];
s=[0  3.4 10.2];
nn=3;
[X,Y]=meshgrid(t,s);
pe=[1E3 0.932e+06 0.967e+06 1.003e+06 1.04e+06 1.106e+06 1.236e+06 1.383e+06 1.55e+06 1.739e+06 1.956e+06 2.206e+06 2.496e+06 2.837e+06 3.245e+06 3.744e+06 4.382e+06 11.022e+06 20.571e+06 31.723e+06 44.337e+06 58.33e+06 388e6 3e10;
1E4 2.08e+06 2.377e+06 2.734e+06 3.172e+06 10.133e+06 21.201e+06 33.945e+06 48.23e+06 63.983e+06 81.166e+06 99.77e+06 119.807e+06 141.312e+06 164.332e+06 188.935e+06 215.199e+06 243.219e+06 273.102e+06 304.967e+06 338.951e+06 375.204e+06 1189e6 9e10;
nn*1E4 nn*2.08e+06 nn*2.377e+06 nn*2.734e+06 nn*3.172e+06 nn*10.133e+06 nn*21.201e+06 nn*33.945e+06 nn*48.23e+06 nn*63.983e+06 nn*81.166e+06 nn*99.77e+06 nn*119.807e+06 nn*141.312e+06 nn*164.332e+06 nn*188.935e+06 nn*215.199e+06 nn*243.219e+06 nn*273.102e+06 nn*304.967e+06 nn*338.951e+06 nn*375.204e+06 nn*1189e6 nn*9e10];
% pe=[1e3 0.932e+06 0.967e+06 1.003e+06 1.04e+06 1.106e+06 1.236e+06 1.383e+06 1.55e+06 1.739e+06 1.956e+06 2.206e+06 2.496e+06 2.837e+06 3.245e+06 3.744e+06 4.382e+06 11.022e+06 20.571e+06 31.723e+06 44.337e+06 58.33e+06 388e6 3e10;
% 1e3 0.932e+06 0.967e+06 1.003e+06 1.04e+06 1.106e+06 1.236e+06 1.383e+06 1.55e+06 1.739e+06 1.956e+06 2.206e+06 2.496e+06 2.837e+06 3.245e+06 3.744e+06 4.382e+06 11.022e+06 20.571e+06 31.723e+06 44.337e+06 58.33e+06 388e6 3e10;
% 1e3 0.932e+06 0.967e+06 1.003e+06 1.04e+06 1.106e+06 1.236e+06 1.383e+06 1.55e+06 1.739e+06 1.956e+06 2.206e+06 2.496e+06 2.837e+06 3.245e+06 3.744e+06 4.382e+06 11.022e+06 20.571e+06 31.723e+06 44.337e+06 58.33e+06 388e6 3e10];
[xq,yq]=deal(Tk,salinity);
varargout{1}= interp2(X,Y,pe,xq,yq,'linear');

end



%---------------relative permeability----------------------%

function varargout = relperm(s, opt, varargin)

kwm = opt.kwm;
s_lmax=1;s_gmax=1-opt.sr(1);
s1=(s(:,1)-opt.sr(1))./(s_lmax-opt.sr(1));
s2=(1-s(:,1)-opt.sr(2))./(s_gmax-opt.sr(2));

s1(s1 < 0) = 0;  s1(s1 > 1) = 1;
s2(s2 < 0) = 0;  s2(s2 > 1) = 1;
varargout{1}    = [ kwm(1) * s1.^2 , kwm(2) * s2.^2];
if nargout > 1
    dkr11=zeros(size(s1,1),1);
    dkr22=zeros(size(s1,1),1);
    ins1=(s1>=0&s1<=1);
    ins2=(s2>=0&s2<=1);
    den1=s_lmax-opt.sr(1);
    den2=s_gmax-opt.sr(2);
    dkr11(ins1)=kwm(1)./den1.*2.*s1(ins1);
    dkr22(ins2)=-kwm(2)./den2.*2.*s2(ins2);
    dkr11(~ins1)=0;
    dkr22(~ins2)=0;
    varargout{2} = [ dkr11,    dkr22];
end




end



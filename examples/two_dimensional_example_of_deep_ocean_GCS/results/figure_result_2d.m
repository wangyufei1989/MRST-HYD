
clear;%clear the workspace
%close all
%close all % close all figures
%% install the mrst_co2
%add the path of the folder of mrst_co2;for windows system, you can use addpath('\..\..'); or simply find the 'startup.m' file in the main folder and run it. 
addpath('../../') 
run startup.m


nx=20;ny=1;nz=40;
dims= [nx ny nz];
% the domain sizes in the (x,y,z) directions are (distance, thickness, depth);
distance=1000;thickness=60; depth=400;
domain=[distance thickness  depth]; 
distance=800;
% generate cartesian grid  system; the dimension can be three;
Grid= computeGeometry(cartGrid(dims,domain));
% change the caresian grid system to radial system
Grid=orth2radial(Grid);
X=Grid.cells.centroids(1:nx);
Y=linspace(0,depth,nz)+3500;
xa=X;
ya=Y;
[X,Y]=meshgrid(X,Y);
%%
s1=0;s2=0.4;
h1=0;h2=0.05;
t1=273.15;t2=276.15+14;
den1=-15;den2=0;
k1=-0.1;k2=0.3;
c1=0;c2=0.07;
m1=0.5;m2=0.54;
%%
figure;
set(gcf,'Position',[20 20 700 750])



%% saturation
dX=0.2;dY=0.108;
dx=0.22;
dy=0.13;
xp=0.05;
yp=0.05;
p1=[xp yp+dy*5+0.03 dX dY];
p2=[xp yp+dy*4 dX dY];

p3=[xp yp+dy*3 dX dY];
p4=[xp yp+dy*2 dX dY];
p5=[xp yp+dy dX dY];
p6=[xp yp dX dY];
%%
p7=[xp+dx yp+dy*5+0.03 dX dY];
p8=[xp+dx yp+dy*4 dX dY];

p9=[xp+dx yp+dy*3 dX dY];
p10=[xp+dx yp+dy*2 dX dY];
p11=[xp+dx yp+dy dX dY];
p12=[xp+dx yp dX dY];

%%
%%
p13=[xp+2*dx yp+dy*5+0.03 dX dY];
p14=[xp+2*dx yp+dy*4 dX dY];

p15=[xp+2*dx yp+dy*3 dX dY];
p16=[xp+2*dx yp+dy*2 dX dY];
p17=[xp+dx*2 yp+dy dX dY];
p18=[xp+dx*2 yp dX dY];
%%
p19=[xp+3*dx yp+dy*5+0.03 dX dY];
p20=[xp+3*dx yp+dy*4 dX dY];

p21=[xp+3*dx yp+dy*3 dX dY];
p22=[xp+3*dx yp+dy*2 dX dY];
p23=[xp+dx*3 yp+dy dX dY];
p24=[xp+dx*3 yp dX dY];

%%
load('s_2dsst5')
subplot('position',p1);

S=reshape(x.s(:,2),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
xlabel('r [m]')
ylabel('z [m]')

caxis([s1 s2])
set (gca,'Ydir','reverse') 
title('t=5 [year]')
%% hydrate volume fraction
%subplot(6,4,5)
subplot('position',p2);
S=reshape(0.3-x.poro(:,1),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([h1 h2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])


%% temperature
%subplot(6,4,9)
subplot('position',p3);
S=reshape(x.Tk(:,1),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([t1 t2])
set (gca,'Ydir','reverse')
xticks([])
yticks([])
%% reaction rate
%subplot(6,4,13)
subplot('position',p4);
S=reshape(((x.r_k(:,1) )) ,nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([k1 k2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
%% co2 concetration
%subplot(6,4,17)
subplot('position',p5);
S=reshape(x.frac(:,2),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([c1 c2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
%% molality of NaCl
%subplot(6,4,21)
subplot('position',p6);
S=reshape(x.m_NaCl(:,1),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([m1 m2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])

%%
%%
load('s_2dst10')
%% saturation
%subplot(6,4,2)
subplot('position',p7);
S=reshape(x.s(:,2),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([s1 s2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
title('t=10 [year]')
%% hydrate volume fraction
%subplot(6,4,6)
subplot('position',p8);
S=reshape(0.3-x.poro(:,1),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([h1 h2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])


%% temperature
%subplot(6,4,10)
subplot('position',p9);
S=reshape(x.Tk(:,1),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])

caxis([t1 t2])
set (gca,'Ydir','reverse')
xticks([])
yticks([])
%% reaction rate
subplot('position',p10);

S=reshape(((x.r_k(:,1) )) ,nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([k1 k2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
%% co2 concetration
subplot('position',p11);
S=reshape(x.frac(:,2),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([c1 c2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
%% molality of NaCl
subplot('position',p12);
S=reshape(x.m_NaCl(:,1),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([m1 m2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
%%
%%
load('s_2dst15')
%% saturation
subplot('position',p13);
S=reshape(x.s(:,2),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([s1 s2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
title('t=15 [year]')
%% hydrate volume fraction
subplot('position',p14);
S=reshape(0.3-x.poro(:,1),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([h1 h2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])

%% temperature
subplot('position',p15);
S=reshape(x.Tk(:,1),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([t1 t2])
set (gca,'Ydir','reverse')
xticks([])
yticks([])
%% reaction rate
subplot('position',p16);
S=reshape(((x.r_k(:,1) )) ,nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([k1 k2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
%% co2 concetration
subplot('position',p17);
S=reshape(x.frac(:,2),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([c1 c2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
%% molality of NaCl
subplot('position',p18);
S=reshape(x.m_NaCl(:,1),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
caxis([m1 m2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
%%
%%
%%
load('s_2dst20')
%% saturation
subplot('position',p19);
S=reshape(x.s(:,2),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])

colorbar
caxis([s1 s2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
title('t=20 [year]')
hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '     \it S_g \rm[-]';
set(colorTitleHandle ,'String',titleString,'fontsize',10);
%% hydrate volume fraction
subplot('position',p20);
S=reshape(0.3-x.poro(:,1),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
colorbar
caxis([h1 h2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '     \phi_{HYD}[-]';
set(colorTitleHandle ,'String',titleString,'fontsize',10);

%% temperature
subplot('position',p21);
S=reshape(x.Tk(:,1),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])

colorbar
caxis([t1 t2])
set (gca,'Ydir','reverse')
xticks([])
yticks([])
hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '    \it T\rm [K]';
set(colorTitleHandle ,'String',titleString,'fontsize',10);
%% reaction rate
subplot('position',p22);
S=reshape(((x.r_k(:,1) )) ,nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
colorbar
caxis([k1 k2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '    \it r_k \rm[mol/s]';
set(colorTitleHandle ,'String',titleString,'fontsize',10);
%% co2 concetration
subplot('position',p23);
S=reshape(x.frac(:,2),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
colorbar
caxis([c1 c2])
set (gca,'Ydir','reverse')
xticks([])
yticks([])
hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '    \it X_l^C \rm[-]';
set(colorTitleHandle ,'String',titleString,'fontsize',10);
%% molality of NaCl
subplot('position',p24);
S=reshape(x.m_NaCl(:,1),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
colorbar
caxis([m1 m2])
set (gca,'Ydir','reverse') 
xticks([])
yticks([])
hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = ' \it m_{\rmNaCl}\rm[molal]';
set(colorTitleHandle ,'String',titleString,'fontsize',10);
colormap jet

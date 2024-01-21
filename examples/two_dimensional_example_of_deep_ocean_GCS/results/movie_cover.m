
clear;%clear the workspace
close all % close all figures
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
Y=linspace(0,depth,nz);
[X,Y]=meshgrid(X,Y);

figure
set(gcf,'color','white')
set(gcf,'position',[553,318,1288,638])
axes('Position',[0.0000003105590062112,0.11,0.995341614906832,0.815])
ax = gca;
ax.NextPlot = 'replaceChildren';
N=10;
A(1:N)=struct('cdata',[],'colormap',[]);




for i=0:20
    savei=num2str(i);
    savex=strcat('s_2dss','t',savei,'.mat');
      subplot(2,1,2)
    if i>0
    load(savex)
    S=reshape(x.s(:,2)+randn(20*40,1)*1e-10,nx,nz);
    end
    if i==0
         S=reshape(randn(20*40,1)*1e-10,nx,nz); 
    end
    contourf(X,Y,S','edgecolor','none')
    %pbaspect([3.5,1,1])
    colormap jet
    AAA= strcat('Saturation of CO_2(g) at t=',savei,'year');
    title(AAA,'fontsize',20)
    colorbar
    set (gca,'Ydir','reverse')
    grid on
    ylabel('Depth below seafloor [m]','fontsize',15)
   xlabel('Radius [m]','fontsize',20)
   
    xlim([ 0 800])
    caxis([0 0.8])
    hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '     S_g[-]';
set(colorTitleHandle ,'String',titleString,'fontsize',20);
    %%
    subplot(2,1,1)
    if i>0
    S=reshape(0.3-x.poro(:,1)+randn(20*40,1)*1e-10,nx,nz);
    end
     if i==0
         S=reshape(randn(20*40,1)*1e-10,nx,nz); 
    end
    contourf(X,Y,S','edgecolor','none')
    %pbaspect([3.5,1,1])
    colormap jet
    AAA= strcat('HYD voulme fraction at t=',savei,'year');
    title(AAA,'fontsize',20)
    colorbar
    set (gca,'Ydir','reverse')
    grid on
    ylabel('Depth below seafloor [m]','fontsize',15)
   xlabel('Radius [m]','fontsize',20)
    xlim([ 0 800])
    caxis([0 0.1])
    hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '     \phi_{HYD}[-]';
set(colorTitleHandle ,'String',titleString,'fontsize',20);
      drawnow
    A(i+1)=getframe(gcf);
    
end

writerObj=VideoWriter('movie_for_sat_hyd.avi');
writerObj.Quality=100;
writerObj.FrameRate=1;
open(writerObj);
writeVideo(writerObj,A);
close(writerObj);




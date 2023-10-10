function dt=Updatedt_kinetic(state,dt)






    dpm=100000*50;
dsm=1e-2*5;

state.ds=abs(state.ds);
state.dp=abs(state.dp);
a=1;

dt=dt*min([5,0.5*abs(dpm/max(max(state.dp(:,1))))*a,0.5*abs(dsm/max(state.ds(:,1)))]);




dt=min([dt, state.CFL*20,40e3]);



end


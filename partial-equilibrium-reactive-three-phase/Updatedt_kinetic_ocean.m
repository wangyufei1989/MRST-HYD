function [dt,state]=Updatedt_kinetic_ocean(state,dt)

if dt==0.1
    state.dtmax=state.dtmax*0.6;
    state.dtmax=max(state.dtmax,1e4);
end

if dt==state.dtmax
    state.dtmax=state.dtmax*1.01;
end




    dpm=100000*20;
dsm=1e-2*5;

state.ds=abs(state.ds);
state.dp=abs(state.dp);
a=1;

dt=dt*min([5,1*abs(dpm/max(max(state.dp(:,1))))*a,1*abs(dsm/max(state.ds(:,1)))]);




dt=min([dt, state.CFL*0.8,state.dtmax]);



end

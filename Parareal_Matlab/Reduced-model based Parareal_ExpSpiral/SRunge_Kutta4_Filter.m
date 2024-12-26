function u=SRunge_Kutta4_Filter(f,t0,t1,u0,n)
[t,u]=Runge_Kutta4_Filter(f,[t0 t1],u0,n);
u=u(end,:);
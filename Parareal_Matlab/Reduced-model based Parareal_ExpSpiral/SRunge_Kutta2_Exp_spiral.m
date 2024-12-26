function u=SRunge_Kutta2_Exp_spiral(f,t0,t1,u0,n)
[t,u]=Runge_Kutta2_Exp_spiral(f,[t0 t1],u0,n);
u=u(end,:);
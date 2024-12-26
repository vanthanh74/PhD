function u=SRunge_Kutta4_Exp_spiral(f,t0,t1,u0,n)
% global c
[t,u]=Runge_Kutta4_Exp_spiral(f,[t0 t1],u0,n);
% u_old = u;
% u(1) = cos(c*(t1-t0))*u_old(1)  - sin(c*(t1-t0))*u_old(2);
% u(2) = sin(c*(t1-t0))*u_old(1)  + cos(c*(t1-t0))*u_old(2);
u=u(end,:);
function u=SExactSol_Exp_spiral(f,t0,t1,u0,n)
[t,u]=ExactSol_Exp_spiral(f,[t0 t1],u0,n);
u=u(end,:);
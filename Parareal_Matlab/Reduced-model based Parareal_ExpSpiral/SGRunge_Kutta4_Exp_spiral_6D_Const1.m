function u=SGRunge_Kutta4_Exp_spiral_6D_Const1(fG,t0,t1,u_init,n)
global c


u0 = u_init;
u_init = u0;



[t,u]=Runge_Kutta4_Exp_spiral(fG,[t0 t1],u_init,n);
u=u(end,:);
u_old0 = u(:,1:6);
u_old1 = u(:,7:12);

%--------------------------------------------------------------------------
E = [-u_old0(1);u_old0(1)-u_old0(2);u_old0(1)-u_old0(3)];
R = [1 0 0; 0 cos(c*(t1-t0)) sin(c*(t1-t0)); 0 -sin(c*(t1-t0)) cos(c*(t1-t0)) ];
R1 = [0 0 0; 0 sin(c*(t1-t0)) 1-cos(c*(t1-t0)); 0 cos(c*(t1-t0))-1 sin(c*(t1-t0)) ];
u = [[u_old0(1:3)';R*u_old0(4:6)'] + 1/c*([u_old1(1:3)';R*u_old1(4:6)'] +[R1*u_old0(4:6)';R1*E])];
u = u';
   
   
   
   
   
   
   
   
   
   
   
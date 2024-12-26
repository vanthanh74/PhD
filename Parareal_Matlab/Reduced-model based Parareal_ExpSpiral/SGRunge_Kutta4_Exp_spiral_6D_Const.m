function u=SGRunge_Kutta4_Exp_spiral_6D_Const(fG1,t0,t1,u_init,n)
global c


u0 =[u_init(1)*cos(t1-t0) + u_init(4)*sin(t1-t0)...
     u_init(2)...
     u_init(3)...
     -u_init(1)*sin(t1-t0) + u_init(4)*cos(t1-t0)...
     u_init(5)...
     u_init(6)    
     ];
u_old0 = u0; 

u1 = [0 ...
      u_init(1)*sin(t1-t0) - u_init(4)*cos(t1-t0) - u_init(3)*(t1-t0) + u_init(4)...
      -u_init(1)*sin(t1-t0) + u_init(4)*cos(t1-t0) + u_init(2)*(t1-t0) - u_init(4)...
      0 ...
      u_init(1)*cos(t1-t0) - u_init(4)*sin(t1-t0) + u_init(6)*(t1-t0) - u_init(1)...
      -u_init(1)*cos(t1-t0) - u_init(4)*sin(t1-t0) - u_init(5)*(t1-t0) + u_init(1)
];
u_old1 = u1;


%--------------------------------------------------------------------------
E = [-u_old0(1);u_old0(1)-u_old0(2);u_old0(1)-u_old0(3)];
R = [1 0 0; 0 cos(c*(t1-t0)) sin(c*(t1-t0)); 0 -sin(c*(t1-t0)) cos(c*(t1-t0)) ];
R1 = [0 0 0; 0 sin(c*(t1-t0)) 1-cos(c*(t1-t0)); 0 cos(c*(t1-t0))-1 sin(c*(t1-t0)) ];
u = [[u_old0(1:3)';R*u_old0(4:6)'] + 1/c*([u_old1(1:3)';R*u_old1(4:6)'] +[R1*u_old0(4:6)';R1*E])];
u = u';
   
   
   
   
   
   
   
   
   
   
   
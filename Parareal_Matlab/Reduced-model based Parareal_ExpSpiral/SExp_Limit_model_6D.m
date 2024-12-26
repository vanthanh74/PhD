function u=SExp_Limit_model_6D(t0,t1,u0,n)
global c
[t,u]=Exp_Limit_model_6D([t0 t1],u0,n);
%  u(1)  = u0(1)+(1/c)*((t1-t0)*(u0(1)+2*u0(2))+sin(c*(t1-t0))*u0(3)+(1-cos(c*(t1-t0)))*u0(4));
%  u(2) = u0(2)+(1/c)*(-(t1-t0)*(2*u0(1)+u0(2))+(cos(c*(t1-t0))-1)*u0(3)+sin(c*(t1-t0))*u0(4));
%  u(3) = cos(c*(t1-t0))*u0(3)+sin(c*(t1-t0))*u0(4)+(1/c)*(-cos(c*(t1-t0))*(2*(t1-t0)*u0(4))+sin(c*(t1-t0))*(2*(t1-t0)*u0(3))+sin(c*(t1-t0))*(2*u0(1)+u0(2))+(1-cos(c*(t1-t0)))*(u0(1)+2*u0(2)) ) ;
%  u(4) = -sin(c*(t1-t0))*u0(3)+cos(c*(t1-t0))*u0(4)+(1/c)*(sin(c*(t1-t0))*(2*(t1-t0)*u0(4))+cos(c*(t1-t0))*(2*(t1-t0)*u0(3))+(cos(c*(t1-t0))-1)*(2*u0(1)+u0(2))+sin(c*(t1-t0))*(u0(1)+2*u0(2)) );  % Limi(t1-t0) model
 u=u(end,:);                    

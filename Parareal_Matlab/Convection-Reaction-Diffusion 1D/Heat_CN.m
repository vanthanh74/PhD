%% Heat equation del(u)/del(t) = del^2(u)/del(x^2) + f(x,t)
%  Crank-Nicolson


function [U1] = Heat_CN(delta_x,delta_t,f,u0,Nx,Nt)
% u_ex=@(x,t) x*(1-x)^2*exp(-2*t);
% f=@(x,t) (4 - 8*x + 4*x^2 - 2*x^3)*exp(-2*t);
% Nx=20;
    
% delta_x = dX;
% 
% delta_t = dT;
% u0 = ux0_restriction;
% Nx = nx_coarse;
% Nt = nt_coarse;


% number_mesh = 6;
% number_mesh_point=zeros(number_mesh,1);
% norm_l2=zeros(number_mesh,1);
% norm_h1=zeros(number_mesh,1);
% for inumber_mesh=1:number_mesh
%     number_mesh_point(inumber_mesh)=Nx;
% 
%     Nt=5;
%     delta_x=1/Nx;
%     delta_t=delta_x^2;
    
    r=delta_t/delta_x^2;
    T=delta_t*Nt;
    I=eye(Nx+1,Nx+1);
    x=zeros(Nx+1,1);
    t=zeros(Nt+1,1);
    for i=1:Nx+1
        x(i)=(i-1)*delta_x;
    end
    for j=1:Nt+1
        t(j)=(j-1)*delta_t;
    end
    %% Creat matrix A
    A=zeros(Nx+1,Nx+1);
    for i=1:Nx+1
        if (i==1)
            A(i,i)=-2*r;
            A(i,i+1)=r;
        elseif(i==Nx+1)
            A(i,i)=-2*r;
            A(i,i-1)=r;
        else
            A(i,i)=-2*r;
            A(i,i+1)=r;
            A(i,i-1)=r;
        end
    end
    
    G=((I-A/2)^-1);
    B=((I-A/2));
 
    %% Creat matrix U0
%     u0=@(x)x*(1-x)^2;
%     U0=zeros(Nx+1,1);
%     for i=1:Nx
%         U0(i)=u0(i);
%     end
    U0 = u0';
    %% Creat matrix U1
    U1=zeros(Nx+1,1);   
    for n=1:Nt+1
        for i=2:Nx
            U1= G*(B*U0 + delta_t*f(x(i),t(n)));
        end
        U0=U1;
    end

    %% Creat matrix u_exact
%     u_exact=zeros(Nx+1,1);
%     for i=1:Nx+1
%         u_exact(i)=u_ex(x(i),T);
%     end
%     figure
%     plot(x,u_exact,'b',x,U1,'r')
%     legend('Exact Solution','Discrete Solution')
%     xlabel('x');ylabel('value');
%     title(['Comparison between exact and discrete solutions with Nx=',num2str(Nx)]);
    %%  Calculate the error on L^2 

%     norm_l2(inumber_mesh)=0;
%     for i=1:Nx+1
%         norm_l2(inumber_mesh)=norm_l2(inumber_mesh)+(U1(i)-u_exact(i))^2*delta_x;
%     end
%     
%     norm_l2(inumber_mesh)=(norm_l2(inumber_mesh))^(1/2);
%     norm_l2(inumber_mesh);
    %% Calculate the error on H1
% 
%     norm_h1(inumber_mesh)=0;
%     for i=1:Nx
%         norm_h1(inumber_mesh)=norm_h1(inumber_mesh)+(((U1(i+1)-u_exact(i+1))-(U1(i)-u_exact(i)))/delta_x)^2*delta_x;
%     end
%     norm_h1(inumber_mesh)=(norm_h1(inumber_mesh))^(1/2);
    %% Refine mesh (increse mesh point) 
%     Nx=2*Nx;
%     Nt=Nt-5;
% end
% %% Figure for errors respect to number of mesh point
% figure
% plot( log(number_mesh_point), -log(norm_l2), 'red',...
%      log(number_mesh_point), -log(norm_h1), 'magenta',...
%     log(number_mesh_point), 2*log(number_mesh_point)-0.75 ,'black');
%    
% xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
% title('Errors');
% legend('norm_l2','norm_h1','2x','Location','NorthEastOutside'); 
end
    
    
    
    

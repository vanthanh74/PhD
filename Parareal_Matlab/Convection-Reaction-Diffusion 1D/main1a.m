%% solve equation -u_xx(x,y)-u_yy(x,y)=f(x,y) with the Dirichlet boundary condition 
clear all
close all
clc
format short
%% Initial informations
ax=0.0;
bx=1.0;
ay=0.0;
by=1.0;
% cases=2;
u_ex = @(x,y)cos(2*pi*x)*sin(2*pi*y^2);
f = @(x,y) 4*pi*cos(2*pi*x)*(pi*sin(2*pi*y^2)*(1+4*y^2)-cos(2*pi*y^2));
N=4;% number of mesh points of first mesh
number_mesh=1;
number_mesh_point=zeros(number_mesh,1);
norm_max=zeros(number_mesh,1);
norm_l2=zeros(number_mesh,1);
norm_maxh1=zeros(number_mesh,1);
norm_h1=zeros(number_mesh,1);

%% Solve discrite solution and refine mesh
for inumber_mesh=1:number_mesh
    
    h=(bx-ax)/N;
    number_mesh_point(inumber_mesh)=N;
    
    %% Create mesh point    
    x=zeros(N+1,1);
    for i_iter=1:N+1
        x(i_iter)=(i_iter-1)*h;
    end
    y=zeros(N+1,1);
    for i_iter=1:N+1
        y(i_iter)=(i_iter-1)*h;
    end
    
    %% Create matrix A 
    A=zeros((N-1)*(N-1),(N-1)*(N-1));
    B=zeros(N-1,N-1);
    I_h=zeros(N-1,N-1);
    for i=1:N-1
        I_h(i,i)=1;
    end
    for i=1:N-1
        if(i==1)
            B(i,i)=4;
            B(i,i+1)=-1;
        else
            if(i==N-1)
                B(i,i)=4;
                B(i,i-1)=-1;
            else
                B(i,i)=4;
                B(i,i-1)=-1;
                B(i,i+1)=-1;
            end
        end
    end
    for i=1:N-1
        if(i==1)
            A((i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1))=B;
            A((i-1)*(N-1)+1:i*(N-1),i*(N-1)+1:(i+1)*(N-1))=-I_h;
        else
            if(i==N-1)
                A((i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1))=B;
                A((i-1)*(N-1)+1:i*(N-1),(i-2)*(N-1)+1:(i-1)*(N-1))=-I_h;
            else
                A((i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1))=B;
                A((i-1)*(N-1)+1:i*(N-1),i*(N-1)+1:(i+1)*(N-1))=-I_h;
                A((i-1)*(N-1)+1:i*(N-1),(i-2)*(N-1)+1:(i-1)*(N-1))=-I_h;
            end
        end
    end
    A=(1/h^2)*A;
    
    %% Create vector F
    F=zeros((N-1)*(N-1),1);
    for i1=1:N-1
        for i2=1:N-1
            F((i1-1)*(N-1)+i2)=f(i1*h,i2*h);
            if(i1==1)
                if(i2==1)
                    F((i1-1)*(N-1)+i2)=F((i1-1)*(N-1)+i2)+u_ex((i1-1)*h,i2*h)/h^2+u_ex((i1)*h,(i2-1)*h)/h^2;
                else
                    if(i2==N-1)
                        F((i1-1)*(N-1)+i2)=F((i1-1)*(N-1)+i2)+u_ex((i1-1)*h,i2*h)/h^2+u_ex((i1)*h,(i2+1)*h)/h^2;
                    else
                        F((i1-1)*(N-1)+i2)=F((i1-1)*(N-1)+i2)+u_ex((i1-1)*h,i2*h)/h^2;
                    end
                end
            else
                if(i1==N-1)
                    if(i2==1)
                        F((i1-1)*(N-1)+i2)=F((i1-1)*(N-1)+i2)+u_ex((i1+1)*h,i2*h)/h^2+u_ex((i1)*h,(i2-1)*h)/h^2;
                    else
                        if(i2==N-1)
                            F((i1-1)*(N-1)+i2)=F((i1-1)*(N-1)+i2)+u_ex((i1+1)*h,i2*h)/h^2+u_ex((i1)*h,(i2+1)*h)/h^2;
                        else
                            F((i1-1)*(N-1)+i2)=F((i1-1)*(N-1)+i2)+u_ex((i1+1)*h,i2*h)/h^2;
                        end
                    end
                else
                    if(i2==1)
                        F((i1-1)*(N-1)+i2)=F((i1-1)*(N-1)+i2)+u_ex((i1)*h,(i2-1)*h)/h^2;
                    else
                        if(i2==N-1)
                            F((i1-1)*(N-1)+i2)=F((i1-1)*(N-1)+i2)+u_ex((i1)*h,(i2+1)*h)/h^2;
                        end
                    end
                end
            end
        end
    end
    
    %% Solve discrete solution
    u=A\F;
    
    
    %% Get exact solution
    u_exact=zeros((N+1),(N+1));
    for i1=1:N+1
        for i2=1:N+1
            u_exact(i1,i2)=u_ex(x(i1),y(i2));
        end
    end
    
    %% Create discrete solution with boundary
     u_dis=u_exact;
    for i1=1:N-1
        for i2=1:N-1
            u_dis(i1+1,i2+1)=u((i1-1)*(N-1)+i2);
        end
    end
    
    
    %% Figure exact and dicrete solutions
    figure
    subplot(1,2,1)
    surf(x,y,u_dis)
    title('u_dis')
%     zlim([0 1])
    subplot(1,2,2)
    surf(x,y,u_exact)
    title(['u_ex, with N=',num2str(N)])
%     zlim([0 1]);
    

%     n=length(u_dis);
 %% Calculate the error on L^infinity
    norm_max(inumber_mesh)=0.0;
    for i=1:N
        for j=1:N
            if (abs(u_dis(i,j)-u_exact(i,j)) > norm_max(inumber_mesh))
                  norm_max(inumber_mesh)=abs(u_dis(i,j)-u_exact(i,j));
            end
        end
    end    
    norm_max(inumber_mesh)
%%  Calculate the error on L^2 

    norm_l2(inumber_mesh)=0;
    for i=1:N
        for j=1:N
            norm_l2(inumber_mesh)=norm_l2(inumber_mesh)+(u_dis(i,j)-u_exact(i,j))^2*h^2;
        end
    end
    
    norm_l2(inumber_mesh)=(norm_l2(inumber_mesh))^(1/2);
    norm_l2(inumber_mesh)
%% Calculate the error on maxH1    

    norm_maxh1(inumber_mesh)=0;
    for i=1:N-1
        for j=1:N-1
            if (abs(((u_dis(i+1,j)-u_exact(i+1,j))-(u_dis(i,j)-u_exact(i,j)))/h) > norm_maxh1(inumber_mesh) )
              norm_maxh1(inumber_mesh)=abs(((u_dis(i+1,j)-u_exact(i+1,j))-(u_dis(i,j)-u_exact(i,j)))/h);
            end
        end
    end
    norm_maxh1(inumber_mesh)

%% Calculate the error on H1

    norm_h1(inumber_mesh)=0;
    for i=1:N-1
        for j=1:N-1
             norm_h1(inumber_mesh)=norm_h1(inumber_mesh)+(((u_dis(i+1,j)-u_exact(i+1,j))-(u_dis(i,j)-u_exact(i,j)))/h)^2*h^2;
        end
    end
    norm_h1(inumber_mesh)=(norm_h1(inumber_mesh))^(1/2);
    
 
%      legend('Exact Solution','Discrete Solution')
%     xlabel('x');ylabel('value');
%       title(['Comparison between exact and discrete solutions with N=',num2str(N)]);
    % ylim([-6 6])
   % xlim([0.1 0.9])
%% Refine mesh (increse mesh point)    
    N=2*N;

end 
%% Figure Norms 
figure
plot(log(number_mesh_point), -log(norm_max),'blue', log(number_mesh_point),-log(norm_l2), 'red',...
    log(number_mesh_point), -log(norm_maxh1), 'cyan', log(number_mesh_point), -log(norm_h1), 'magenta',...
    log(number_mesh_point), 2*log(number_mesh_point),'black');
    
xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Errors');
legend('norm_max','norm_l2','norm_maxh1','norm_h1','2x','Location','NorthEastOutside');  

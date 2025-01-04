clear;
close all;
clc;
n=10;
c=30;
% Sign_Matrix=ones(n);
global A;
% A=zeros(n);
% sign_vector=[1, 2*(rand(1,n-1)<0.5)-1];
% Sign_Matrix(1,:)=sign_vector;
% A(1,2:n)=rand(1,n-1);
% opt_lambda=c*rand(1,n);
% A(1,1)=sum(A(1,2:n))+opt_lambda(1);
% for i=2:n
%     Sign_Matrix(i,(i+1):n)=-sign_vector(i)*sign_vector((i+1):n);
%     A(i,(i+1):n)=rand(1,n-i);
%     Sign_Matrix(i,1:(i-1))=Sign_Matrix(1:(i-1),i)';
%     A(i,1:(i-1))=A(1:(i-1),i)';
%     A(i,i)=sum(A(i,:))+opt_lambda(i);
% end
% A=A.*Sign_Matrix;
% for i=1:n
%         if sum(abs(A(i,1:(i-1))))+sum(abs(A(i,(i+1):n)))>A(i,i)
%             disp("False")
%         end
%         if abs(opt_lambda(i)-(A(i,i)-(sum(abs(A(i,1:(i-1))))+sum(abs(A(i,(i+1):n))))))>1e-10
%             disp("False")
%         end
% end
% for i=1:n
% for j=(i+1):n
% for k=(j+1):n
% if A(i,j)*A(i,k)*A(j,k)>0
% disp("False")
% end
% end
% end
% end
lambda_min_opt_B_M=1;
while(lambda_min_opt_B_M>1e-4)
A=rand(n);
A=A*A'+c*diag(rand(1,n));
fun = @(Lambda)-trace(diag(Lambda));
Lambda0=zeros(1,n);
lb=zeros(1,n);
ub=[];
AA = [];
b = [];
Aeq = [];
beq = [];
nonlcon=@psdcon;
[opt_lambda_M,opt_trace_M]=fmincon(fun,Lambda0,AA,b,Aeq,beq,lb,ub,nonlcon);
opt_trace_M=-opt_trace_M;
Opt_Lambda=diag(opt_lambda_M);
Opt_B=A-Opt_Lambda;
lambda_min_opt_B_M=min(eig(Opt_B))
end
% % Time Graph
% G=gsp_ring(n,1);
% S=full(G.L);
% [U,D]=eig(S);
% Cov_x=U*A*U';
% Cov_z=U*Opt_Lambda*U';
% Cov_w=U*Opt_B*U';
% mu_x=zeros(1,n);mu_z=zeros(1,n);
% Cov_x_z=[Cov_x,Cov_z;Cov_z, Cov_z];
% x_z=mvnrnd([mu_x,mu_z],Cov_x_z);
% x=x_z(1:n);z=x_z((n+1):2*n);
% estimated_z=Cov_z*inv(Cov_x)*x';
% plot(0:(n-1), x, 'k', 0:(n-1), z, 'b', 0:(n-1), x-z, 'g', 0:(n-1), estimated_z, 'r');
% xlabel("Time");ylabel("Process Values"); legend("x Realization","z Realization","w Realization", "Estimated z Given x");
% title("Stationary Approximation of a Time Process");
% % 2D-Grid Graph
% G = gsp_2dgrid(m);
% Random Sensor Graph
param.distribute=1;
G = gsp_sensor(n,param);
paramplot.show_edges = 1;
paramplot.vertex_size = 200;
paramplot.climits = [-7 7];
% gsp_plot_graph(G,paramplot);
S=full(G.L);
[U,D]=eig(S);
Cov_x=U*A*U';
Cov_z=U*Opt_Lambda*U';
Cov_w=U*Opt_B*U';
mu_x=zeros(1,n);mu_z=zeros(1,n);
Cov_x_z=[Cov_x,Cov_z;Cov_z, Cov_z];
x_z=mvnrnd([mu_x,mu_z],Cov_x_z);
x=x_z(1:n);z=x_z((n+1):2*n);
estimated_z=Cov_z*inv(Cov_x)*x';
figure;
sgtitle("Stationary Approximation of a Process on a General Random Sensor Graph");
subplot(2,2,1);
gsp_plot_signal(G,x,paramplot);
title("x Realization");
subplot(2,2,2);
gsp_plot_signal(G,z,paramplot);
title("z Realization");
subplot(2,2,3);
gsp_plot_signal(G,x-z,paramplot);
title("w Realization");
subplot(2,2,4);
gsp_plot_signal(G,estimated_z,paramplot);
title("Estimated z Given x");
E_x=trace(A)
E_z=trace(Opt_Lambda)
E_w=trace(Opt_B)
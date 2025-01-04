clear;
close all;
clc;
n=100;
c=1000;
Sign_Matrix=ones(n);
global A;
A=zeros(n);
sign_vector=[1, 2*(rand(1,n-1)<0.5)-1];
Sign_Matrix(1,:)=sign_vector;
A(1,2:n)=rand(1,n-1);
opt_lambda=c*rand(1,n);
A(1,1)=sum(A(1,2:n))+opt_lambda(1);
for i=2:n
    Sign_Matrix(i,(i+1):n)=-sign_vector(i)*sign_vector((i+1):n);
    A(i,(i+1):n)=rand(1,n-i);
    Sign_Matrix(i,1:(i-1))=Sign_Matrix(1:(i-1),i)';
    A(i,1:(i-1))=A(1:(i-1),i)';
    A(i,i)=sum(A(i,:))+opt_lambda(i);
end
A=A.*Sign_Matrix;
for i=1:n
        if sum(abs(A(i,1:(i-1))))+sum(abs(A(i,(i+1):n)))>A(i,i)
            disp("False")
        end
        if abs(opt_lambda(i)-(A(i,i)-(sum(abs(A(i,1:(i-1))))+sum(abs(A(i,(i+1):n))))))>1e-10
            disp("False")
        end
end
for i=1:n
for j=(i+1):n
for k=(j+1):n
if A(i,j)*A(i,k)*A(j,k)>0
disp("False")
end
end
end
end
Opt_Lambda=diag(opt_lambda);Opt_B=A-Opt_Lambda;
% Time Graph
G=gsp_ring(n,1);
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
plot(0:(n-1), x, 'k', 0:(n-1), z, 'b', 0:(n-1), x-z, 'g', 0:(n-1), estimated_z, 'r');
xlabel("Time");ylabel("Process Values"); legend("x Realization","z Realization","w Realization", "Estimated z Given x");
title("Stationary Approximation of a Time Process");
E_x=trace(A)
E_z=trace(Opt_Lambda)
E_w=trace(Opt_B)
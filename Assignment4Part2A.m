clear
clc


%Define the initial Parameters
s = 2; 
n = 1; 
x  = sym('x',[n,1]);
syms t real
u = 5;
dtRK4 = 0.1;

%define test funtion
lambda = -2;
f = @(x,t)(lambda*x);


% Define K as a matrix nxs
K=sym('K%d%d',[n,s]);

%Butcher Tableau Definition

A = [1/4 (1/4 - sqrt(3)/6);
    (1/4+sqrt(3)/6) 1/4 ;];      
b = [1/2;1/2;];
c = [(1/2 - sqrt(3)/6) (1/2 + sqrt(3)/6)];

% Create Symbolic Function
X = sym('X%d',[n,1]);
syms t dt real

% IRK4 Scheme Algorithm using Newton iteration
r = sym(zeros(n,s));
for i = 1 : s 
    
    r(:,i) = f(X+dt*K*A(i,:)', t+c(i)*dt)-K(:,i);

end
r=reshape(r,n*s,1);
K2 = reshape(K,n*s,1);
dr = jacobian(r,K2);
matlabFunction(r,dr, 'file', 'rFunction','vars',{t,X,dt,K2});
clear K x dt t


%Initial Conditions Definition

dt          = 0.1; 
tolerance   = 1e-6; 
x0=1;
lambda = -2;
tFinal  = 2;
N  = round(tFinal/dt);
t  = zeros(1,N+1);
xk      = zeros(length(x0),N+1);
xk(:,1) = x0;


%Implicit RK4%

K = zeros(n,s);
for z =1:s 
        K(:,z)=xk(:,1);        
end
for k = 1:N    

    K=reshape(K,n*s,1);
  
    [r,dr] = rFunction(t(k),xk(:,k),dt,K);
    
    while norm(r) > tolerance
        dK = -(dr\r);
        K = K + dK; 
        [r,dr] = rFunction(t(k),xk(:,k),dt,K);
    end
    K = reshape(K,n,s); 
    
    %Record x(k+1) 
     xk(:,k+1) = xk(:,k) + dt*K*b;
     t(k+1) = t(k)+dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Explicit Runge Kutta Method of order 4 %%%

a21= 1/2;
a32 = 1/2;
a43 = 1/2;

b1 = 1/6;
b2 = 1/3;
b3 = 1/3;
b4 = 1/6;

c2 = 1/2;
c3 = 1/2;
c4 = 1;

nRK4 = tFinal / dtRK4;

tRK4 = zeros(nRK4,1);
xRK4 = zeros(nRK4,1);

for j = 2:nRK4+1
    tRK4(j) = tRK4(j-1) + dtRK4;
    lambda_RK4(1) = lambda; 
    xRK4(1) = x0;
    
    K1 = f( xRK4(j-1) );
    K2 = f( xRK4(j-1)+a21*dtRK4*K1 );
    K3 = f( xRK4(j-1)+a32*dtRK4*K2 );
    K4 = f( xRK4(j-1)+a43*dtRK4*K3 );
    xRK4(j) = xRK4(j-1)+dtRK4*((b1*K1)+(b2*K2)+(b3*K3)+(b4*K4));
end

%plot both the Implicit,Explicit and true solution RK4 results
figure (1)
for subfig = 1:n
    title('Implicit and Explicit RK4 Comparison Plot')
    subplot(n,1,subfig);hold on
    plot(t,xk(subfig,:),'r')
    plot(tRK4,xRK4,'marker','.','markersize',10)
    xlabel('t');
    ylabel('x');
    legend('Implicit RK4','Explicit RK4')
end
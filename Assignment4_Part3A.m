clear all
close all
clc


% Parameters for Lagrange Fully Implicit DAE Matlab Function definition
s = 2; 
n = 6; 
M = 1;
g = 9.81;
m = 1;
L=1;

K=sym('K%d%d',[n,1],'real');
x=sym('x%d%d',[n,1],'real');

syms t real
syms z real
 
 f =[x(4:6)-K(1:3);-M*g*[0;0;1]-z*x(1:3)-K(4:6);x(1:3)'*(-M*g*[0;0;1]-z*x(1:3))+x(4:6)'*x(4:6)];

 % Determine Matlab function with four input and output DAEfunc
matlabFunction(f, 'file', 'DAEfunc','vars',{x,z,t,K});
clear x t f z 




% Determine K as matrix nxs and z as matrix mxs using symbolic function
K=sym('K%d%d',[n,s]);
z=sym('z%d%d',[m,s]);

% Butcher Tableau

A = [1/4 (1/4 - sqrt(3)/6);
    (1/4+sqrt(3)/6) 1/4 ;];      
b = [1/2;1/2;];
c = [(1/2 - sqrt(3)/6) (1/2 + sqrt(3)/6)];

% r and dr Matlab Function Definition
X = sym('X%d',[n,1]);
syms t dtIRK real

r = sym(zeros(n+m,s));

for i = 1 : s 
    
    r(:,i) = DAEfunc(X+dtIRK*K*A(i,:)',z(i), t+c(i)*dtIRK,K(:,i));

end
r=reshape(r,(n+m)*s,1);
K=reshape(K,n*s,1);
z=reshape(z,m*s,1);
V = [K;z];
dr = jacobian(r,V);

%Determine Matlab function as having 5 input variables
%(t,X,dtIRK,K,z) and two output r and dr
matlabFunction(r,dr, 'file', 'rFunc','vars',{t,X,dtIRK,K,z});
clear K x dtIRK t z

%Initial Condition
dtIRK    = 0.0001; 
Tol   = 1e-6; 

lambda = -2;
tFinal  = 2;
Nsteps  = round(tFinal/dtIRK);
t  = zeros(1,Nsteps+1);

x0=[1;0;0;0;1;0];
xk      = zeros(length(x0),Nsteps+1);
xk(:,1) = x0;

K = zeros(n,s);
z = zeros(m,s);
C = zeros(1,Nsteps+1);

%Start of the IRK Algorithm
for i =1:s 
        K(:,i)=xk(:,1);        
end


for k = 1:Nsteps    
   
    z=reshape(z,m*s,1);
    K=reshape(K,n*s,1);
    V = [K;z];
     
    [r,dr] = rFunc(t(k),xk(:,k),dtIRK,K,z);
    
    while norm(r) > Tol
        dV = -(dr\r);
        V = V + dV; 
        [r,dr] = rFunc(t(k),xk(:,k),dtIRK,V(1:n*s),V(n*s+1:end));
    end
    
    K=reshape(V(1:s*n,:),n,s);
    z=reshape(V(s*n+1:end,:),m,s); 
     
    %Record the constraint
    C(k)=1/2*(xk(1:3,k)'*xk(1:3,k)-L^2);
    xk(:,k+1) = xk(:,k) + dtIRK*K*b;
    t(k+1) = t(k)+dtIRK;
end



%plot the Dynamics and the Constraints
figure(1)
plot3(xk(1,:), xk(2,:), xk(3,:));
xlabel('x')
ylabel('y')
zlabel('z')

figure(2)
plot(t,C)
xlabel('t (Seconds)')
ylabel('C (DAE Constraint)')
grid;


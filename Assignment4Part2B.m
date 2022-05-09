clear

%Define Initial Condition
s = 2; 
n = 2;
x  = sym('x',[n,1]);
syms t real
u = 5;

% Define Van Der Pol
f = [x(2);
    u*(1-x(1)^2)*x(2) - x(1)];


matlabFunction(f, 'file', 'VanDerPol','vars',{x,t});
clear x t f


% Define K as a matrix size nxs
K=sym('K%d%d',[n,s]);

%State Butcher Tableau

A = [1/4 (1/4 - sqrt(3)/6);
    (1/4+sqrt(3)/6) 1/4 ;];      
b = [1/2;1/2;];
c = [(1/2 - sqrt(3)/6) (1/2 + sqrt(3)/6)];

% Define residual function and find its jacobian derivative
X = sym('X%d',[n,1]);
syms t dt real
%K = reshape(K,n*s,1);
r = sym(zeros(n,s));
for i = 1 : s 
    
    r(:,i) = VanDerPol(X+dt*K*A(i,:)', t+c(i)*dt)-K(:,i);

end

r=reshape(r,n*s,1);
K2 = reshape(K,n*s,1);
dr = jacobian(r,K2);

% Name the function with 4 input and 2 output
matlabFunction(r,dr, 'file', 'rFunction2','vars',{t,X,dt,K2});
clear K x dt t

% Initial Condition definition
dt          = 1e-2; 
tolerance   = 1e-6; 

x0=[2;0];
tFinal  = 25;
N  = round(tFinal/dt);
t  = zeros(1,N+1);

xk      = zeros(length(x0),N+1);
xk(:,1) = x0;


% IRK Scheme Algorithm using Newton iteration
K = zeros(n,s);

%Guess for K value
for z =1:s 
        K(:,z)=xk(:,1);        
end



for k = 1:N    

    K=reshape(K,n*s,1);
    
    
    [r,dr] = rFunction2(t(k),xk(:,k),dt,K);
    
    
    while norm(r) > tolerance
        dK = -(dr\r);
        K = K + dK; 
        [r,dr] = rFunction2(t(k),xk(:,k),dt,K);
    end
    K = reshape(K,n,s); 
    
    %Record x(k+1) value 
     xk(:,k+1) = xk(:,k) + dt*K*b;
     t(k+1) = t(k)+dt;
end


%Plot the Implicit Runge-Kutta order 4
figure (1)
for subfig = 1:n
    title('VanDerPol Solution using Implicit Runge-Kutta 4')
    subplot(n,1,subfig);hold on
    plot(t,xk(subfig,:),'ro')
    xlabel('Time t')
    ylabel('Solution (x,y)')
end




clc;clear;
close all

%% Defining the function:
lambda = -2;
f = @(x)(lambda*x);

%% Simulation parameters:
tFinal  = 2;
dtEuler = 0.1;
dtRK2   = 0.1;
dtRK4 = 0.1;
x0 = 1;

%% The exact Solution:
tExact = 0:0.1:tFinal;
xExact = x0*exp(lambda*tExact);

%% Explicit Euler Method:

nEuler = tFinal / dtEuler;

t_Euler = zeros(nEuler,1);
xEuler = zeros(nEuler,1);
 
 
for j = 2:nEuler+1
    t_Euler(j) = t_Euler(j-1) + dtEuler;
    
    xEuler(1) = x0;
    xEuler(j) = xEuler(j-1) + dtEuler*f(xEuler(j-1));
end

%% Runge Kutta Method of order 2:
a21= 1/2;
b2 = 1;
b1 = 0;
c2 = 1/2;

nRK2 = tFinal / dtRK2;

tRK2 = zeros(nRK2,1);
xRK2 = zeros(nRK2,1);


for j = 2:nRK2+1
    tRK2(j) = tRK2(j-1) + dtRK2;
% RK2 step
        a21= 1/2;
        b1 = 0;
        b2 = 1;
        c2 = 1/2;
        
        xRK2(1) = x0;
        
       K1_RK2 = f( xRK2(j-1));
       K2_RK2 = f( xRK2(j-1)+ a21*dtRK2*K1_RK2);
       xRK2(j) = xRK2(j-1) + dtRK2*a21*K1_RK2 + dtRK2 * c2 * K2_RK2 ;
end

%% Runge Kutta Method of order 4
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




%% Compare the results:
figure(1)
plot(tExact,xExact,'marker','.','markersize',20)
hold on
plot(t_Euler,xEuler,'marker','.','markersize',10)
plot(tRK2,xRK2,'marker','.','markersize',10)
plot(tRK4,xRK4,'marker','.','markersize',10)
xlabel('t');
ylabel('x');
legend('Exact','Euler Method','RK2','RK4')


%% Integration for various time steps:
N_table = floor(logspace(1,3,50));

for index = 1:length(N_table)
    
    N = N_table(index);
    dt(index) = tFinal / N_table(index);
    t = zeros(N,1);
    xEuler = zeros(N,1);
    xRK2 = zeros(N,1);
    xRK4 = zeros(N,1);
    
   
    
    for j = 2:N+1
       
        % Euler step
        b1 = 1;
        
        xEuler(1) = x0;
       
         
        xEuler(j) = xEuler(j-1) + dt(index)* b1 * f(xEuler(j-1));
        
    end
    
    
   for j = 2:N+1
       
        % RK4 step
        a21= 1/2;
        a32 = 1/2;
        a43 = 1;

        b1 = 1/6;
        b2 = 1/3;
        b3 = 1/3;
        b4 = 1/6;
        
        c2 = 1/2;
        c3 = 1/2;
        c4 = 1;
        
        xRK4(1) = x0;
        
    K1 = f( xRK4(j-1));
    K2 = f( xRK4(j-1)+ c2*dt(index)*K1);
    K3 = f( xRK4(j-1)+ c3*dt(index)*K2);
    K4 = f( xRK4(j-1)+ c4*dt(index)*K3);
    xRK4(j) = xRK4(j-1)+ dt(index)*((b1*K1) + (b2*K2) + (b3*K3) + (b4*K4));
   end
   
    for j = 2:N+1
        
        % RK2 step
        a21= 1/2;
        b2 = 1;
        c2 = 1/2;
        
        xRK2(1) = x0;
        
       K1_RRK2 = f( xRK2(j-1));
       K2_RRK2 = f( xRK2(j-1)+ a21*dt(index)*K1_RRK2);
       xRK2(j) = xRK2(j-1) + dt(index) * a21 * K2_RRK2 + dt(index) * c2 * K2_RRK2 ;
       
    end
    
    
    % Compute the global error
    err(1,index) = norm(xEuler(end)-xExact(end),inf);
    err(2,index) = norm(xRK2(end)-xExact(end),inf);
    err(3,index) = norm(xRK4(end)-xExact(end),inf);
   
end

figure(2);clf; 
loglog(dt,err,'marker','.','markersize',15,'linestyle','none')
grid on
set(gca,'XDir','reverse')
xlabel('Delta t')
ylabel('Global Error (|| X_N-X(T) ||)')
legend('Euler', 'RK2','RK4')

%%  Van der Pol Oscillator Task 3

% Assume initial x = 2, y = 0
[t,y] = ode45(@Vanderpol_order_5,(0:0.01:25),[2; 1]);

% Plot the solutions for  and  against t.
figure (3)
YY = [y(:,1),y(:,2)];
stem(t,YY)
title('Van der Pol (u = 5) with ODE45');
xlabel('Time t');
ylabel('Solution (x,y)');
% xline(5.28,11.09,16.89,22.7'-',{'Asymptote 1= 5.28');
legend('x','y')
% 


% Runge-Kutta iteration 



dtRK4_v = 0.01 ;
tFinal_v = 25;

nRK4_v = tFinal_v / dtRK4_v;

tRK4_v = zeros(nRK4_v,1);
xRK4_v = zeros(2,nRK4_v);
xRK4_v(1) = 2;
% 
for j = 2:nRK4_v+1
    
        a21= 1/2;
        a32 = 1/2;
        a43 = 1;

        b1 = 1/6;
        b2 = 1/3;
        b3 = 1/3;
        b4 = 1/6;
        
        c2 = 1/2;
        c3 = 1/2;
        c4 = 1;
    
         tRK4_v(j) = tRK4_v(j-1) + dtRK4_v;
         
    
         
    K1_1 = Vanderpol_order_5( tRK4_v(j-1), xRK4_v(:,j-1));   
    K2_1 = Vanderpol_order_5( tRK4_v(j-1)+a21*dtRK4_v, xRK4_v(:,j-1)+c2*dtRK4_v*K2); 
    K3_1 =Vanderpol_order_5( tRK4_v(j-1)+a32*dtRK4_v, xRK4_v(:,j-1)+c3*dtRK4_v*K3);  
    K4_1 = Vanderpol_order_5( tRK4_v(j-1)+a43*dtRK4_v, xRK4_v(:,j-1)+c4*dtRK4_v*K4);   
   
    xRK4_v(:,j) = xRK4_v(:,j-1)+dtRK4_v.*((b1.*K1_1)+(b2.*K2_1)+(b3.*K3_1)+(b4.*K4_1));
    xRK4_vv =xRK4_v.';
end

% % Plot the solutions for  and  against t.
figure (4)
X = [xRK4_vv(:,1),xRK4_vv(:,2)];
stem(tRK4_v,X)
title('Van der Pol (u = 5) with Runge-Kutta 4th order');
xlabel('Time t');
xlim([0 25]);
ylabel('Solution (x,y)');
% xline(5.28,11.09,16.89,22.7'-',{'Asymptote 1= 5.28');
legend('x','y','Location','northwest')
% 

% Function Definition of Van Der Pol Oscillator wherre Y = [x y]
function dydt = Vanderpol_order_5(t,Y)

dydt = [Y(2); 5*(1-Y(1)^2)*Y(2)-Y(1)];

end


% 
% 
% 


% 
% 

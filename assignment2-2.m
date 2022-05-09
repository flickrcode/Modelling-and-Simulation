load input.mat 
load output.mat

%  1) Model definition 

% In this exercise we will work with the 3 candidate ARX defined as 
% y(t)+a1*y(t-1)+a2*y(t-2) = b0*u(t)+e(t)
% y(t)+a1*y(t-1)+a2*y(t-2) = b0*u(t)+b1*u(t-1)+e(t)
% y(t)+a1*y(t-1)+a2*y(t-2)+a3*y(t-3) = b1*u(t-1)+e(t)



%%
% 2. Identification using Least Squares formula
% First of all lets split the data in estimation and validation sets
% (half and half)
N = length(y); 

% number of data
uest = u(1:N/2);
yest = y(1:N/2);
uval = u(N/2+1:end);
yval = y(N/2+1:end);


% Build the matrix PHI:
PHI = zeros(3,N/2);           % The matrix has 3 rows (3 parameters), N/2 columns (time instants)
% PHI = zeros(4,N/2);         % Switch to this PHI for the 2nd and 3rd Candidate

% Switch for each respective candidate
PHI(:,1) = [ uest(1) ; 0; 0 ]; %candidate 1
% PHI(:,1) = [ uest(1) ; uest(2); 0; 0]; %candidate 2
% PHI(:,1)   = [ uest(1) ; 0; 0; 0 ];  %candidate 3

for i=N/2:-1:3 % change the index to i=N/2:-1:4 for the 3rd ARX
    
    % Switch for each respective candidate
    PHI(:,i) = [uest(i) ; yest(i-1); yest(i-2)] ;                   %candidate 1
%     PHI(:,i) = [uest(i) ; uest(i-1); yest(i-1); yest(i-2)];      %candidate 2
%       PHI(:,i) = [uest(i-1) ; yest(i-1); yest(i-2); yest(i-3)];  %candidate 3
end

% Write the expression in code as well, so we get our estimate!
      th = (PHI*PHI')\PHI*yest; 

% th should be be a vector with two elements, estimated value for a and
% estimate value for b. Pick this two elements separetely:

% Switch for each respective candidate
%First ARX Parameters
bhat_0 = th(1)
ahat_1 = -th(2)
ahat_2 = -th(3)

%Second ARX Parameters
% bhat_0 = th(1)
% bhat_1 = th(2)
% ahat_1 = -th(3)
% ahat_2 = -th(4)

%Third ARX Parameters
% bhat_1 = th(1)
% ahat_1 = -th(2)
% ahat_2 = -th(3)
% ahat_3 = -th(4)


%% 3. ARX Prediction & Simulation

N = length(y);
NN = N/2;

% Redefine the data variables so it will be easy to switch from
% validation to estimation sets.
yn = yval;
un = uval;
% yn = yest; 
% un = uest;



% Write your code for 1-step-ahead predictor
ypred = zeros(NN,1); % the vector where we will store the predicted output

for i=NN:-1:3 % For the 3rd ARX model, we switch the index to i=NN:-1:4
    
    %Switch between the 1st ARX, 2nd ARX and 3rd ARX for the predictor
    %candidate 1
    ypred(i) = -ahat_1*yn(i-1)-ahat_2*yn(i-2)+ bhat_0*un(i); 
    %candidate 2
%   ypred(i) = -ahat_1*yn(i-1)-ahat_2*yn(i-2) + bhat_0*un(i) +bhat_1*un(i-1);
    %candidate 3
%   ypred(i)= -ahat_1*yn(i-1)-ahat_2*yn(i-2)-ahat_3*yn(i-3) + bhat_1*un(i);
end


% Write the code for simulation
ysim = zeros(NN,1); % the vector where we will store the simulated output

for i=NN:-1:3 %For the 3rd ARX model, we switch the index to i=NN:-1:4
    
     %Switch between the 1st ARX, 2nd ARX and 3rd ARX for the simulation
     %candidate 1
    ysim(i) = -ahat_1*ysim(i-1)-ahat_2*ysim(i-2)+ bhat_0*un(i);
    %candidate 2
%   ysim(i) = -ahat_1*ysim(i-1)-ahat_2*ysim(i-2) + bhat_0*un(i) +bhat_1*un(i-1);
    %candidate 3
%   ysim(i)= -ahat_1*ysim(i-1)-ahat_2*ysim(i-2)-ahat_3*ysim(i-3) + bhat_1*un(i-1);
end

% compare with real data and compute RMSE
predERROR = yn-ypred;
predRMSE = rms(predERROR)
simERROR = yn-ysim;
simRMSE  = rms(simERROR)


% 4. Covariance Analysis of the ARX models
%Gaussian White Noise
noise_std = 0.01;
covariance = inv(PHI*PHI.')*(1/NN)*(noise_std^2)
% plot Predicted DATA vs MODEL prediction and 
figure (1)

subplot(2,1,1)
plot(yn)
hold on
plot(ypred)
legend(' DATA','Model prediction')
title('Output')
xlabel('Samples')
ylabel('output')
subplot(2,1,2)
plot(predERROR)
legend('prediction error')
xlabel('Samples')
ylabel('error')

figure (2)

subplot(2,1,1)
plot(yn)
hold on
plot(ysim)
legend('DATA','Model simulation')
title('Output')
xlabel('Samples')
ylabel('output')
subplot(2,1,2)
plot(simERROR)
legend('Simulation error')
xlabel('Samples')
ylabel('error')


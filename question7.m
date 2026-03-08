clc
clear
close
% NOTE: last section takes a couple of minutes to run

%% 7.1 finding estimator for alpha and sigma refining dt
% Parameters
rng(30026); %Random seed
alpha = 1; %Exact alpha
sigma = 1; %exact sigma
X0 = 1; %Initial condition
T = 1e3; %Final time
h = 2^(-10); %Fine timestep for "exact" solution
nSteps = T/h; %Number of steps for exact solution

% Generating exact solution single trajectory
Y = zeros(nSteps+1, 1);
Y(1) = X0;
dW = sqrt(h)*randn(1, nSteps); %brownian increments
for ii = 1:nSteps
    Y(ii+1) = Y(ii) - alpha * Y(ii) * h + sqrt(2*sigma) * dW(ii); %Euler-Maruyama method
end

% Computing alpha estimator as analytically found in Q6
dt = 2.^-(0:7); %sampling rate values
nSample = T./dt; % Number of samples for each sampling rate
alpha_tilde = zeros(length(dt),1);
for ii = 1:length(nSample)
    indeces = 1 : nSteps/nSample(ii) : nSteps+1;
    X_tilde = Y(indeces);
    alpha_tilde(ii) = -1/dt(ii) * log( sum(X_tilde(1:end-1) .* X_tilde(2:end)) / sum(X_tilde(1:end-1).^2) ); %expression found in Q6.4
end

%Plotting results
figure
semilogx(dt, alpha_tilde, '-or', dt, alpha*ones(length(dt),1), '-b', 'LineWidth',2)
title('Alpha Estimator')
legend('alpha estimator', 'alpha real')
xlabel('dt')
ylabel('alpha')

%% 7.1bis finding estimator for alpha and sigma pushing T to infinite
%The result in Q6 shows that the converges is reached when N goes to
%infinity, independently of the sampling rate. Hence we fix a sampling rate
%and try to vary the values of T in order to have bigger values of N
rng(30026); %Random seed
alpha = 1; %Exact alpha
sigma = 1; %exact sigma
X0 = 1; %Initial condition
T = 1e5; %Final time
h = 2^(-10); %Fine timestep for "exact" solution
nSteps = T/h; %Number of steps for exact solution

% Generating exact solution single trajectory
Y = zeros(nSteps+1, 1);
Y(1) = X0;
dW = sqrt(h)*randn(1, nSteps); %brownian increments
for ii = 1:nSteps
    Y(ii+1) = Y(ii) - alpha * Y(ii) * h + sqrt(2*sigma) * dW(ii); %Euler-Maruyama method
end

% Computing alpha estimator as analytically found in Q6
dt = 2^-(4); %fixed sampling rate value
Tstop_vec = [1e2, 1e3, 1e4, 1e5]; %increasing final time of the samplings
alpha_tilde = zeros(length(Tstop_vec),1);
for ii = 1:length(Tstop_vec)
    Tstop = Tstop_vec(ii);
    indeces = 1 : dt/h : nSteps * Tstop/T + 1; 
    X_tilde = Y(indeces);
    alpha_tilde(ii) = -1/dt * log( sum(X_tilde(1:end-1) .* X_tilde(2:end)) / sum(X_tilde(1:end-1).^2) ); %expression found in Q6.4
end

%Plotting results
figure
semilogx(Tstop_vec, alpha_tilde, '-or', Tstop_vec, alpha*ones(length(Tstop_vec),1), '-b', 'LineWidth',2)
title('Alpha Estimator')
legend('alpha estimator', 'alpha real')
xlabel('T')
ylabel('alpha')

%% 7.2 verifying that the estimators satisfy the central limit theorem
rng(30026);
M = 1e4; %number of simulations
alpha = 1;
sigma = 1;
h = 2^(-9); %number of timesteps
X0 = 1; %initial position
T = 1e3; %final time
nSteps = T/h; 
dt = 1; %sample rate
nSample = T/dt; %number of samples
alpha_tilde = zeros(M,1); %initializing estimator for alpha
%computing estimator for alpha
for m = 1:M
    Y = zeros(nSteps+1, 1);
    Y(1) = X0;
    dW = sqrt(h)*randn(1, nSteps);
    for ii = 1:nSteps
        Y(ii+1) = Y(ii) - alpha * Y(ii) * h + sqrt(2*sigma) * dW(ii);
    end
    indeces = 1 : nSteps/nSample : nSteps+1; %indeces for the sample
    X_tilde = Y(indeces);
    alpha_tilde(m) = -1/dt * log(sum( X_tilde(1:end-1) .* X_tilde(2:end) ) / sum( X_tilde(1:end-1).^2) );
end

%plotting distribution of the estimator and normal distribution
var = (exp(2*dt*alpha) - 1)/dt^2; %variance for normal distribution
figure
histogram(sqrt(nSample)*(alpha_tilde-alpha), 'Normalization', 'pdf'); 
x = linspace(min(-4*sqrt(var)), max(4*sqrt(var)), 100);
normal_pdf = normpdf(x, 0, sqrt(var));
hold on;
plot(x, normal_pdf, 'r', 'LineWidth', 2);
hold off;
title('Estimator central limit theorem')
legend('Histogram', 'Normal distribution')

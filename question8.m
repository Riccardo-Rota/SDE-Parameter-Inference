clc
clear
close

%% 8 Numerical Results for eigenvalues-based estimators for alpha and sigma

% Explicit nonlinear system
g1 = @(a,s,X,delta) sum( X(1:end-1).^2 .* (X(2:end) - exp(-a*delta).X(1:end-1) + X(2:end).^2 - s/a - exp(-2*a*delta)(X(1:end-1).^2 - s/a)) );
g2 = @(a,s,X,delta) sum( X(1:end-1) .* (X(2:end) - exp(-a*delta).X(1:end-1) + X(2:end).^2 - s/a - exp(-2*a*delta)(X(1:end-1).^2 - s/a)) );

% Parameters
rng(30026); %Random seed
alpha = 1; %Exact value of alpha
sigma = 1; % Exact value of sigma
X0 = 1; %Initial condition
T = 1e3; %Final time
h = 2^(-10); %Fine timestep for "exact" solution
nSteps = T/h; % Number of steps for exact solution

% Generating exact solution single trajectory
Y = zeros(nSteps+1, 1);
Y(1) = X0;
dW = sqrt(h)*randn(1, nSteps); %Generating brownian increments
for ii = 1:nSteps
    Y(ii+1) = Y(ii) - alpha * Y(ii) * h + sqrt(2*sigma) * dW(ii); %Euler-Maruyama method
end
dt = 1; %Fixed sampling rate
nSample = T/dt; %Number of evaluated points on the trajectory
indeces = 1 : nSteps/nSample : nSteps+1; %Generating the grid of sample points
X_tilde = Y(indeces); %Generating syntethic observations
asol=zeros(nSample,2); %Initializing vector of estimator pairs
% Computing estimators varying number of observed values
for ii = 2:nSample
    G = @(est) [g1(est(1),est(2),X_tilde(1:ii), dt), g2(est(1),est(2),X_tilde(1:ii), dt)];
    init_guess = [0.5, 0.5];
    asol(ii, :) = fsolve(G, init_guess); %Solving the system
end

% Plotting results for alpha
figure
plot(1:nSample, asol(:, 1));
hold on
yline(1, '-r', '\alpha = 1', 'LineWidth', 1.5)
title('\alpha estimators')
xlabel('number of available observations')

% Plotting results for sigma
figure
plot(1:nSample, asol(:, 2));
hold on
yline(1, '-r', '\sigma = 1', 'LineWidth', 1.5)
title('\sigma estimators')
xlabel('number of available observations')
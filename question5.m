clc
clear
close

%% computing classic estimators for alpha and sigma
rng(30026);

% Experimental setting
alpha = 1;
sigma = 1;
X0 = 1; %initial value
T = 1e3; %final time
h = 2^(-10); %timestep for the "exact" trajectory
nSteps = T/h;
tVec = 0:h:T;

% Generating "exact" single trajectory
Y = zeros(nSteps+1, 1);
Y(1) = X0;
dW = sqrt(h)*randn(1, nSteps); %random brownian increments
for ii = 1:nSteps
    Y(ii+1) = Y(ii) -alpha * Y(ii) * h + sqrt(2*sigma) * dW(ii); %Milstein method
end

% Computing alpha and sigma as definied in question4, with different
% sampling rates
dt = 2.^-(0:7); %vector of progressively finer sampling rates
N = T./dt;
%initializing vectors to be populated in the for loop
alpha_hat = zeros(length(dt),1); 
sigma_hat = zeros(length(dt),1);
%computing the estimators for different values of dt
for ii = 1:length(N)
    indeces = 1 : nSteps/N(ii) : nSteps+1; %indeces of the sampled points
    X_tilde = Y(indeces); %sampled points
    alpha_hat(ii) = -1/dt(ii) * sum(X_tilde(1:end-1) .* (X_tilde(2:end)- X_tilde(1:end-1))) / sum(X_tilde(1:end-1).^2);
    sigma_hat(ii) = 1/(2*N(ii)*dt(ii)) * sum( (X_tilde(2:end)- X_tilde(1:end-1)).^2 );
end

%plot of the estimator for alpha
figure(1)
semilogx(dt, alpha_hat, '-or', dt, alpha*ones(length(dt),1), '-b', 'LineWidth',2)
title('Alpha Estimator')
legend('alpha estimator', 'alpha real')
xlabel('\Delta')
ylabel('alpha')

%plot of the estimator for sigma
figure(2)
semilogx(dt, sigma_hat, '-or', dt, sigma*ones(length(dt),1), '-b', 'LineWidth',2)
title('Sigma Estimator')
legend('sigma estimator', 'sigma real')
xlabel('\Delta')
ylabel('sigma')

% Plot of the trend of the error on the estimate of alpha
err_alpha = alpha*ones(length(dt)) - alpha_hat;
figure(3)
semilogx(dt, dt/2, '-b', dt, err_alpha, '-or', 'LineWidth',2)
title('Error of alpha estimator')
legend('\Delta/2', 'error')
xlabel('\Delta')
ylabel('error')

% Plot of the trend of the error on the estimate of sigma
err_sigma = sigma*ones(length(dt)) - sigma_hat;
figure(4)
semilogx(dt, dt/2, '-b', dt, err_sigma, '-or', 'LineWidth',2)
title('Error of sigma estimator')
legend('\Delta/2', 'error')
xlabel('\Delta')
ylabel('error')
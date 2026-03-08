clc
clear
close

%% 2 plotting distribution of final time solution

% Parameters of the problems
rng(30026); %Random seed
T = 1e3; %Final timestep
M = 1e4; %Number of realizations
h = 2^(-9); %Discretization step
nSteps = T/h; %Number of timesteps
alpha = 1; % Exact drif coefficient
sigma = 1; %Exact diffusion coefficient
X0 = 1; %Inizial condition
% Generating final time for all the trajectories
Y_final = zeros(M, 1);
for m = 1:M
    Y = X0;
    dW = sqrt(h)*randn(1, nSteps); %Generating brownian increments
    for ii = 1:nSteps
        Y = Y - alpha * Y * h + sqrt(2*sigma) * dW(ii); %Euler-Maruyama scheme
    end
    Y_final(m) = Y; %Saving final timestep
end
%Plotting the distribution of Y_final and mu_inf
figure
histogram(Y_final, 'Normalization', 'pdf');
x = linspace(min(-4), max(4), 100);
normal_pdf = normpdf(x, 0, sigma/alpha);
hold on;
plot(x, normal_pdf, 'r', 'LineWidth', 2);
hold off;
title('Final time distribution of numerical solution')
legend('Y final', '\rho_{\infty}')

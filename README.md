# Parameter Inference for Stochastic Differential Equations (SDEs)

## Overview
This repository investigates the problem of inferring the parameters of a Stochastic Differential Equation (SDE) given discrete-time observations of a single continuous trajectory. Specifically, the project focuses on deriving, implementing, and analyzing estimators for the drift ($\alpha$) and diffusion ($\sigma$) coefficients of the Ornstein-Uhlenbeck (OU) process. 

The OU process is widely used in physics and biology, and its governing SDE is defined over $t\in[0,T]$ as:

$$dX_{t}=-\alpha X_{t}dt+\sqrt{2\sigma}dW_{t}$$

## Mathematical Framework & Methodology 
The research compares two primary approaches to parameter estimation:

### 1. Classic Estimators
Standard discrete-time approximations of the SDE are heavily dependent on the sampling rate $\Delta$. Mathematical analysis and numerical simulations reveal that these classic estimators are inherently biased. Their expected values only converge to the true parameters in the limit as the sampling rate approaches zero ($\Delta\rightarrow0$).

### 2. Eigenvalue-Based (Spectral) Estimators
In real-world applications, data is often provided at a fixed sampling rate, making the limit $\Delta\rightarrow0$ impossible to achieve. To solve this, this repository implements a more sophisticated estimating strategy based on the eigenvalues and eigenfunctions of the SDE's generator operator, $-\mathcal{L}_{a}$. 

Unlike classic estimators, these eigenvalue-based estimators are strongly consistent; they converge to the true parameters almost surely as the number of observations $N\rightarrow\infty$, completely independent of the fixed sampling rate $\Delta$.

## Repository Structure & MATLAB Scripts
The codebase is implemented in MATLAB. To reproduce the experiments, run the numbered scripts sequentially:

* **`01_stationary_distribution.m`**: Simulates the OU process using the Euler-Maruyama method. It verifies the ergodicity of the solution by comparing the final-time numerical distribution of $10^4$ realizations against the theoretical stationary density $\rho_{\infty}$.
* **`02_classic_estimators.m`**: Computes the classic estimators for $\alpha$ and $\sigma$ across progressively finer sampling rates $\Delta=2^{-i}$. It outputs convergence plots demonstrating that the estimation error decreases linearly with $\Delta$.
* **`03_eigenvalue_estimator_alpha.m`**: Evaluates the eigenvalue-based estimator for the drift coefficient $\alpha$ assuming $\sigma$ is known. It empirically proves the estimator's independence from $\Delta$ and verifies that the estimator satisfies the Central Limit Theorem (CLT).
* **`04_eigenvalue_estimators_joint.m`**: Solves the complex case where both $\alpha$ and $\sigma$ are unknown. It constructs a two-dimensional non-linear system based on the first two eigenfunctions and solves it iteratively using `fsolve` over an increasing number of observations.

## Read the Full Report
For full mathematical proofs—including martingale property checks, stationary covariance calculations, and Taylor expansions of the estimator errors—please see the attached `Report_SDE_Inference.pdf`.

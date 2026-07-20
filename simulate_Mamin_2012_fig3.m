
% Simulation of Surface Spin Correlations (Mamin et al. 2012)
clear; clc; close all;


gamma = 2 * pi * 28.024e9; % Gyromagnetic ratio for NV (rad/(s*T))
tau = 8e-6;                % Total echo time (s)
tD = linspace(0, tau, 100); % Pulse position within echo (s)

% Contrast parameters from Section VI
alpha0 = 1;                % Normalized bright counts
alpha_tau = 0.91;          % Measured alpha(tau)/alpha0
beta_tau = 0.70;           % Measured beta(tau)/alpha0 beta0 is the dark counts


% Format: [T1 (s), Brms (Tesla)]% Specific cases from Fig 3 in Mamin paper
cases = [
    20e-6, 500e-9;  % Case 1: 50 us, 370 nT
    50e-6, 500e-9;  % Case 2: 20 us, 425 nT
    1000e-6,  500e-9;  % Case 3: 6 us,  500 nT
    10000000e-6,  500e-9   % Case 4: 2 us,  700 nT
];
labels = {'T1=20us, 500nT', 'T1=50us, 500nT', 'T1=100us, 500nT', 'T1=10 s, 500nT'};


figure(1); hold on;
for i = 1:size(cases, 1)
    T1 = cases(i, 1);
    Brms = cases(i, 2);
    
    % Equation (6): Dimensionless factor f(tau, tD, T1)
    %  |tau/2 - tD| is the term used in the paper for pulse offset
    offset = abs(tau/2 - tD);
    
    f = (2*T1/tau^2) * (tau - 5*T1 + 4*T1*exp(-offset/T1) ...
        - 2*T1*exp(-(tau/2 + offset)/T1) + 2*T1*exp(-(tau/2 - offset)/T1) ...
        + T1*exp(-tau/T1));
    
    % Equation (7): Fluorescence Signal
    S_avg = (alpha_tau + beta_tau) / 2;
    S_diff = (alpha_tau - beta_tau) / 2;
    S = S_avg + S_diff * exp(-(gamma^2 * Brms^2 * tau^2 .* f) / 2);
    
    plot(tD * 1e6, S, 'LineWidth', 1.5);
end

xlabel('t_D (\mus)');
ylabel('Signal (a.u.)');
%title('Simulation of Fig 3: Effect of Pulse Position within Hahn Echo');
legend(labels);

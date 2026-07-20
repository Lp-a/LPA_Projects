
% Load experimental data
[B_exp, spc_exp] = eprload('filename.DTA');
% Normalize experimental data (makes fitting easier)
spc_exp = spc_exp / max(abs(spc_exp));
% Convert B_exp from Gauss (G) to milliTesla (mT) since EasySpin works in mT
B_exp = B_exp / 10;


%%%% SIMULATION %%%%%%

% Sys structure holds the spin system parameters
Sys.S = 1/2; % Electron spin (S=1/2 for nitroxide)
Sys.Nucs = '14N'; % Nuclear spin (I=1 for 14N)

% g-tensor (MHz/mT conversion factor is 2.8025 GHz/T, or 28.025 MHz/mT)
% Typical g-tensor for TEMPOL (axial symmetry, from literature):
Sys.g = [2.0076, 2.0076, 2.0024]; % [g_x, g_y, g_z]
% Isotropic g-value for quick calculation:
g_iso = sum(Sys.g) / 3; % 2.005867
% Conversion factor: mu_B / h approx 28.025 MHz/mT
CONVERSION_FACTOR = g_iso * 28.025; 


% A-tensor (Hyperfine Coupling, usually in MHz)
% Typical A-tensor for TEMPOL (axial symmetry, from literature, in MHz):
% A_par ~ 90-100 MHz, A_perp ~ 15-20 MHz
A_par = 105.0; % Parallel component (MHz) 96
A_perp = 17.0; % Perpendicular component (MHz) 17
Sys.A = [A_perp, A_perp, A_par]; % [A_x, A_y, A_z]

% Line Broadening Parameters
% lwpp (linewidth peak-to-peak) is the width of the individual spectral lines.
% This represents the residual, unresolved broadening (homogenous/inhomogenous).
Sys.lwpp = 0.1; % Gaussian broadening (G) use 0.15 (0.05) if broadened from instrumental factors, unresolved proton hyperfine, etc.

%% 2. Define Experimental Conditions
% from experiments
Exp.mwFreq = 9.844267; % Microwave frequency (GHz) - Typical X-band frequency 339 mT
Exp.Range = [340 365]; % Field sweep range (mT) - Adjust to center around the spectrum
Exp.nPoints = length(B_exp); % Number of points in the simulated spectrum

% To fit for tau_c - the rotational diffusion constant (D).
% Will fit log(D) as it spans many orders of magnitude.
% D = 1 / (6 * tau_c)
Initial_tau_c = 6e-9; % Initial guess for tau_c (seconds)
Initial_D = 1 / (6 * Initial_tau_c);

% Set up the System structure for chili with the initial guess
SysFit = Sys;
SysFit.Diff = Initial_D; % Initial D
SysFit.Liouvilletype = 'fastest'; % Required solver for nitroxides

% Setup the Fit structure
Vary.Diff = 0.5; % Allow Diff to vary by 50% around the initial guess
Vary.lwpp = 0.1; % Allow lwpp to vary by 0.1 mT
Vary.g = [0 0 0]; % Fixing g-tensor
Vary.A = [0 0 0]; % Fixing A-tensor

% Define the range limits for the variables (min, max)
% D is 1/(6*tau_c). tau_c from 1e-10 to 1e-7 s corresponds to D from 1.6e8 to 1.6e5 rad^2/s
% D_min (tau_c=1e-7s) = 1.6e5; D_max (tau_c=1e-10s) = 1.6e8
% We will fit log10(D)
Par0.logDiff = log10(Initial_D); % Start fitting in log space
ParMax.logDiff = 8.2; % Max log10(D) (tau_c ~ 1e-9 s)
ParMin.logDiff = 5.0; % Min log10(D) 5 (tau_c ~ 1e-6 s)

% Set the initial residual broadening value
Par0.lwpp = Sys.lwpp;
ParMax.lwpp = 0.5;
ParMin.lwpp = 0.05;

% Options for the fit function
FitOpt.Method = 'leastsq'; % Or 'leastsq' 'nelder-mead'(Levenberg-Marquardt)
FitOpt.Verbosity = 1;
FitOpt.OutlierThr = 0.2; % Exclude points far from the current simulation

%% 4. Perform Manual Fit by Iterating over tau_c
% Define the range of tau_c values to test (Logarithmically spaced is best)
tau_c_min = 1e-9; % 0.5 ns (Faster limit)
tau_c_max = 15e-9; % 10 ns (Slower limit)
num_tests = 50; % Number of points to test

tau_c_test_list = logspace(log10(tau_c_min), log10(tau_c_max), num_tests);
rmsd_list = zeros(1, num_tests); % To store the goodness-of-fit (RMSD or chi-squared)

% Use a fixed, reasonable lwpp for the initial fit
fixed_lwpp = 0.15; % You can also try fitting lwpp manually later

SysFit = Sys; % Start with the rigid-limit system
SysFit.lwpp = fixed_lwpp; 
SysFit.Liouvilletype = 'fastest'; % Required solver for nitroxides

% Find the experimental B range indices to match the simulation's Exp.Range
% NOTE: Ensure B_exp and Exp.Range are compatible!
idx = (B_exp >= Exp.Range(1)) & (B_exp <= Exp.Range(2));
B_exp_match = B_exp(idx);
spc_exp_match = spc_exp(idx);

fprintf('\n--- Starting Manual Tau_c Scan (%d points) ---\n', num_tests);
for i = 1:num_tests
    tau_c = tau_c_test_list(i);
    
    % Calculate Rotational Diffusion Constant (D)
    D = 1 / (6 * tau_c);
    
    % Simulate the spectrum using CHILI
    SysFit.Diff = D; 
    [~, spec_sim] = chili(SysFit, Exp);
    
    % Truncate the simulation spectrum to match the experimental B-field points
    spec_sim_match = spec_sim(idx); 
    
    % Normalize the simulated spectrum
    spec_sim_match = spec_sim_match / max(abs(spec_sim_match));
    
    % Calculate RMSD (Root-Mean-Square Deviation) for goodness of fit
    residuals = spc_exp_match - spec_sim_match;
    
    % CRITICAL FIX: Use 'all' to ensure a scalar mean output
    rmsd = sqrt(mean(residuals.^2, 'all')); 
    
    rmsd_list(i) = rmsd; % This line will now assign a scalar value
    
    if mod(i, 10) == 0
        fprintf('  Tested tau_c = %.2e s. RMSD = %.4f\n', tau_c, rmsd);
    end
end
fprintf('--- Scan Complete ---\n');

% Find the tau_c that gives the minimum RMSD
[min_rmsd, min_idx] = min(rmsd_list);
tau_c_fit = tau_c_test_list(min_idx);

% Re-simulate with the best-fit tau_c for plotting
D_fit = 1 / (6 * tau_c_fit);
SysFit.Diff = D_fit;
[B_fit, SpecFit] = chili(SysFit, Exp);
SpecFit = SpecFit / max(abs(SpecFit));

%% 5. Plot Results

% Plotting RMSD vs. tau_c
figure(3);
semilogx(tau_c_test_list, rmsd_list, 'ko-', 'LineWidth', 1.5);
hold on;
plot(tau_c_fit, min_rmsd, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Rotational Correlation Time, $\tau_c$ (s)', 'Interpreter', 'latex');
ylabel('RMSD (Goodness-of-Fit)');
title('Manual Fit: RMSD vs. \tau_c Scan');
grid on;

% Plotting the Best Fit vs. Experiment
figure(4);
plot(B_exp, spc_exp, 'k-', 'LineWidth', 1.5);
hold on;
% Need to ensure the simulation is plotted on the experimental B_exp axis
plot(B_fit(idx), SpecFit(idx), 'r--', 'LineWidth', 1.5); 
xlabel('Magnetic Field (mT)');
ylabel('EPR Signal (1st Derivative)');
title('Best Manual Fit vs. Experimental Data');
legend('Experimental Data', sprintf('Fitted Simulation (\\tau_c = %.2e s)', tau_c_fit), 'Location', 'best');
grid on;

fprintf('\n--- Manual Fit Results ---\n');
fprintf('Best Fitted Rotational Correlation Time (tau_c): %.2e seconds\n', tau_c_fit);
fprintf('Minimum RMSD: %.4f\n', min_rmsd);
fprintf('-------------------\n');

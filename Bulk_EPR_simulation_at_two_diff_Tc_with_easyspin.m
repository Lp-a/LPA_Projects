clear;
close all;

% Sys structure holds the spin system parameters
Sys.S = 1/2; % Electron spin (S=1/2 for nitroxide)
Sys.Nucs = '14N'; % Nuclear spin (I=1 for 14N)

% g-tensor (MHz/mT conversion factor is 28.025 MHz/mT)
% Typical g-tensor for TEMPOL (axial symmetry, from the JPC B paper omim cl values):
Sys.g = [2.0094, 2.0061, 2.0019]; % [g_x, g_y, g_z]
% Isotropic g-value for quick calculation:
% g_iso = (g_x + g_y + g_z) / 3 = 2.0058
g_iso = sum(Sys.g) / 3; % 2.005867
% Conversion factor: mu_B / h approx 28.025 MHz/mT
CONVERSION_FACTOR = g_iso * 28.025; 

% A-tensor (Hyperfine Coupling, usually in MHz)
% Typical A-tensor for TEMPOL (from the JPC B paper omim cl values, in MHz):
% A_par ~ 90-100 MHz, A_perp ~ 15-20 MHz
% A(MHz) = A(G) * gamma_e/1e6
% gamma_e (Hz/G) ≈ 2.802495e6 Hz/G
A_zz = 95.845329; % Parallel component (MHz)
A_xx = 19.3372155; % Perpendicular component (MHz)
A_yy = 15.693972; % Perpendicular component (MHz)
Sys.A = [A_xx, A_yy, A_zz]; 

% Line Broadening Parameters
% lwpp (linewidth peak-to-peak) is the width of the individual spectral lines.
% This represents the residual, unresolved broadening (homogenous/inhomogenous).
Sys.lwpp = 0.15; % Gaussian broadening (G)


% Exp structure holds the experimental settings
Exp.mwFreq = 9.5; % Microwave frequency (GHz) - Typical X-band frequency 339 mT
Exp.Range = [332 345]; % Field sweep range (mT) - Adjust to center around the spectrum
Exp.nPoints = 1024; % Number of points in the simulated spectrum


% Liq structure holds the rotational motion parameters for fast motion usic
% garlic and Chili for slow motion from 1-100 ns Tau_c 
% Used the isotropic model for this simple simulation.
% Tau_c is related to viscosity and temperature (Stokes-Einstein-Debye theory).

% Fast Motion (low viscosity, uses Redfield Theory)
% tau_c = 5e-11 s (50 ps)
LiqFast.tcorr = 130e-12; % Rotational correlation time (seconds)
LiqFast.Diff = 1/(6*LiqFast.tcorr); % Rotational diffusion constant (rad^2/s) - for isotropic motion

% Intermediate Motion (High viscosity, typical RTIL at room temp, uses Stochastic Liouville Theory)
% The paper suggests TEMPOL in RTILs might be in the intermediate regime (1e-10 to 1e-8 s)
tau_c_rtil = 40e-9; % 3 ns
D_rtil = 1/(6*tau_c_rtil);
OptRTIL_chili.tcorr = tau_c_rtil; % Rotational correlation time (seconds)



% Fast Motion 
[B1, spec1] = garlic(Sys, Exp, LiqFast); % 'garlic' is for liquid-phase simulation
spec1 = spec1 / max(spec1); % Normalize

% Intermediate Motion (RTIL) 
SysRTIL = Sys; % Copy the full ANISOTROPIC Sys
SysRTIL.Diff = D_rtil; % Set the isotropic diffusion constant (D, in rad^2/s)
SysRTIL.Liouvilletype = 'fastest'; % Fast SLE solution for nitroxides

[B2, spec2] = chili(SysRTIL, Exp);
spec2 = spec2 / max(spec2); % Normalize

% Plotting
figure(1);
plot(B1, spec1, 'b-', 'LineWidth', 1.5);
hold on;
plot(B2, spec2, 'r--', 'LineWidth', 1.5);
xlabel('Magnetic Field (mT)');
ylabel('EPR Signal (1st Derivative)');
%title('Simulated TEMPOL EPR Spectrum');
legend(sprintf('Fast Motion (\\tau_c = %.1e s)', LiqFast.tcorr), ...
       sprintf('RTIL Motion (\\tau_c = %.1e s)', OptRTIL_chili.tcorr), ...
       'Location', 'best');
grid on;
set(gca, 'FontName', 'Helvetica', 'FontSize', 12);


tau_c_start = 0.2e-9; 
tau_c_end = 10e-9;   
num_points = 8;

tau_c_list = logspace(log10(tau_c_start), log10(tau_c_end), num_points);% use this to create evenly spaced across different oreders of mag
%tau_c_list = linspace(tau_c_start, tau_c_end, num_points);
linewidth_list = zeros(1, num_points); 
 
ExpLinewidth = Exp;
ExpLinewidth.Range = [334 337]; % Adjust this if the peak to peak linewidth estimation fails

for i = 1:num_points
    tau_c = tau_c_list(i);
    
   
    if tau_c < 1e-9 
        % Use GARLIC
        Liq.tcorr = tau_c;
        [B, spec] = garlic(Sys, ExpLinewidth, Liq);
    else 
        % Use CHILI
        SysLoop = Sys; 
        D = 1/(6*tau_c);
        SysLoop.Diff = D; 
        SysLoop.Liouvilletype = 'fastest';
        [B, spec] = chili(SysLoop, ExpLinewidth);
    end
    

    % Find the index of the global max (first peak of the first derivative)
    [~, max_idx_all] = max(spec);
     max_idx = max_idx_all(1);% Ensure only the first index is used, even if min returns multiple
    
    % Find the index of the minimum immediately following the max
    % Max is restricted to the data points after the minimum.
     [~, min_rel_idx_all] = min(spec(max_idx+1:end));
    
 
    if isempty(min_rel_idx_all)
        % If no peak is found after the max, skip
        linewidth_list(i) = NaN;
        fprintf('Warning: Could not find minimum after maximum for tau_c = %.1e s. Assigning NaN.\n', tau_c);
        continue; % Skip to the next iteration
    end
    
    min_rel_idx = min_rel_idx_all(1);
    
    % Convert the relative index back to the absolute index in the 'B' vector
    min_idx = max_idx + min_rel_idx;
    
    % Calculate the linewidth (distance between min and max)
    linewidth = B(min_idx) - B(max_idx);
    
    % Assign the linewidth value
    linewidth_list(i) = linewidth;
    linewidth_freq_list(i) = linewidth_list(i) * CONVERSION_FACTOR;
end

% Plotting the Linewidth vs. Log(tau_c)
figure(2);
semilogx(tau_c_list, linewidth_freq_list, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
xlabel('Rotational Correlation Time (s)', 'Interpreter', 'latex');
ylabel('Linewidth (MHz)');
%title('EPR Linewidth Dependence on Rotational Correlation Time');
grid on;
set(gca, 'FontName', 'Helvetica', 'FontSize', 12);

disp('Linewidth vs. Tau_c analysis complete.');

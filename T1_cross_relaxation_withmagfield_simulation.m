% Constants
mu0 = 4*pi*1e-7;         % vacuum permeability (T*m/A)
muB = 9.274e-24;         % Bohr magneton (J/T)
g = 2;                   % g-factor of free electron
gamma_e = 2*pi*2.8e6;    % rad/s/G
B0 = 300;                % magnetic field in Gauss
D = 2.87e9*2*pi;         % zero-field splitting (rad/s)
w0 = D - gamma_e*B0;     % NV transition lower frequency (rad/s)

% NV depth / hemisphere radius
r_NV = 10e-9;            % 10 nm
Nspins = 1e6;           % number of spins

% Spin density
V_hemi = (2/3)*pi*r_NV^3;
n = Nspins / V_hemi;     % spins/m^3

% Magnetic field variance
B2 = (8*pi/9) * n * (mu0/(4*pi))^2 * (g*muB)^2 / (r_NV^3);

% Dipolar and rotational correlation times
tau_dip = 1 / (mu0*(g*muB)^2*n/1.054e-34);  % simplified estimate
tau_rot = 1e-8;  % example, set for molecular rotation at ~10^8 Hz
tau_c = 1 / (1/tau_dip + 1/tau_rot);

% Noise spectral density and T1
S_B = (B2 * tau_c) / (1 + (w0*tau_c)^2);

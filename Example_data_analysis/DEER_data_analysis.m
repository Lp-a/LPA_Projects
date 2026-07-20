
%path = '';

nits = 20;%548
npts = 26;%26
t_tot = 520; % ms
t = linspace(20,t_tot,npts);

window_radius = 20;
window = [window_radius window_radius];

sig_mat = zeros(npts,nits);
ref_mat = zeros(npts,nits);

for jj = 1:nits
    
    currfile = ['Iteration_' num2str(jj) '.mat'];
    load([path,currfile])

    sig_mat(:,jj) = savedData.signal;
    ref_mat(:,jj) = savedData.reference;

end

contrast_mat = (ref_mat-sig_mat)./(ref_mat);
mean_contrast_cum = contrast_mat(:, end);%contrast_mat(:, end)
%mean_contrast_cum = mean(contrast_mat, 2);

% Define the model function
decaying_sine = @(b, t) b(1) * exp(-t / b(2)) .* cos(2*pi*b(3)*t + b(4)) + b(5);

% Initial guesses for parameters: [A, tau, f, phi, offset]
b0 = [0.1, 200, 0.01, 0, 0];

% Fit using nonlinear least squares
opts = optimset('Display', 'off');
[b_fit, ~, residuals, exitflag] = lsqcurvefit(decaying_sine, b0, t(:), mean_contrast_cum(:), [], [], opts);

% Generate fitted curve
t_fit = linspace(min(t), max(t), 500);
fit_curve = decaying_sine(b_fit, t_fit);

% Plot data and fit
figure;
plot(t, mean_contrast_cum, 'DisplayName', 'Data');
hold on;
plot(t_fit, fit_curve, '-', 'LineWidth', 2, 'DisplayName', 'Fit');
xlabel('Time (ms)');
ylabel('Mean Contrast');
legend('Location', 'best');
title('Decaying Sine Fit to Mean Contrast');
grid on;


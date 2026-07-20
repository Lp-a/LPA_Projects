close all
clear

%path  = 'C:\Users\lpa3\Box\Backlund lab shared\data\Lakshmy\2026\March\3rd\single in MEA new magnet position\2026_March3_xy2_NV2_DEER_T1_sushkov_xyy_pi_xyy_1pulse\';
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\Lakshmy\2026\March\10th\2026_March10_xy8_NV1_DEER_T1_sushkov_1_pulse/';
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\Lakshmy\2026\March\10th\xy4_NV1_DEER_T1_xyy_pi_xyy_xyy_xyy_1_pulse\';
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\Lakshmy\2026\March\9th\singles in MEA\2026_March9_xy5_NV1_DEER_T1_xyy_pi_xyy__xyy_xyy_1_pulse\';
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\Lakshmy\2026\April\4th\2026_April4_xy1_NV1_40mW_DEER_T1_sushkov_with_noise_0dbm_2p5Vpp\';

%path  = 'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\2026_April7_xy2_NV1_DEER_T1_sushkov_noise_neg5dbm_2p5Vpp\';
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\2026_April6_xy4_NV1_DEER_T1_sushkov_noise\';
%path =  'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\2026_April3_xy3_NV1_DEER_T1_sushkov_with_noise_neg5dbm\';
path = 'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\2026_April3_xy3_NV1_DEER_T1_sushkov_with_noise\';
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\2026_April5_xy1_NV1_DEER_T1_sushkov_with_noise\';

%nits  = 37 ;
%nits  = 55;
%nits = 45;
%nits  = 100;
%nits  = 35;

%nits = 100;%200
%nits = 20;%100
%nits = 200;%275
nits = 40; %67
%nits = 35; %55

npts = 26;

t_tot = 15;%us
window_radius = 4;
window = [window_radius window_radius];

t = linspace(0.01,t_tot,npts);


sig_mat = zeros(npts,nits);
ref_mat = zeros(npts,nits);


for jj = 1:nits
    
    currfile = ['Iteration_' num2str(jj) '.mat'];
    load([path,currfile])

    sig_mat(:,jj) = savedData.signal;
    ref_mat(:,jj) = savedData.reference;

end

sig_mat_separateruns = zeros(size(sig_mat));
ref_mat_separateruns = zeros(size(sig_mat));

sig_mat_separateruns(:,1) = sig_mat(:,1);
ref_mat_separateruns(:,1) = ref_mat(:,1);

for jj = 2:nits
    
    sig_mat_separateruns(:,jj) = jj*sig_mat(:,jj) - (jj-1)*sig_mat(:,jj-1);
    ref_mat_separateruns(:,jj) = jj*ref_mat(:,jj) - (jj-1)*ref_mat(:,jj-1);


end


sig_mat_separateruns = cumsum(sig_mat_separateruns,2)./repmat(1:nits,[npts,1]);
ref_mat_separateruns = cumsum(ref_mat_separateruns,2)./repmat(1:nits,[npts,1]);


contrast_cum = (ref_mat(:,nits)-sig_mat(:,nits))./(ref_mat(:,nits));
c_min_cum = min(contrast_cum);
c_max_cum = max(contrast_cum);
%norm_contrast_cum = (contrast_cum - c_min_cum) / (c_max_cum - c_min_cum);
norm_factor = mean(contrast_cum(1));
norm_contrast_cum  = contrast_cum / norm_factor;

contrast_mat = (ref_mat_separateruns-sig_mat_separateruns)./(ref_mat_separateruns);
contrast_mat_norm = (contrast_mat - c_min_cum) / (c_max_cum - c_min_cum);
%contrast_mat_norm = contrast_mat / norm_factor;

std_err_mean = std(contrast_mat,[],2)/sqrt(size(contrast_mat,2));
norm_std_err = std(contrast_mat_norm, [], 2) / sqrt(size(contrast_mat_norm, 2));
norm_std_dev = std(contrast_mat_norm, [], 2);
P95 = tinv(0.95, size(contrast_mat,2)-1);
CI95 = mean(contrast_mat,2) + std_err_mean*[-1 1]*P95;

ft = fittype('a + b*exp(-c*x)'); % 'exp(-x/tau)'  'a + b*exp(-c*x)'
options = fitoptions(ft);

%options.StartPoint = 5; % Initial guess for tau us

options.StartPoint = [0, 1, 1]; 
options.Lower = [-0.5, 0, 0.02]; 
options.Upper = [0.5, 1.5, 50];


% Perform fit
[curve, gof] = fit(t', norm_contrast_cum, ft, options);


%alpha = 1 - 0.50 = 0.50
%ci50 = paramci(pd, 'Alpha', 0.50);

ci50 = confint(curve,0.5);

c_val   = curve.c;
c_lower = ci50(1,3); 
c_upper = ci50(2,3);

T1_val   = 1 / c_val;
T1_high_bound = 1 / c_lower;
T1_low_bound  = 1 / c_upper;

err_high = T1_high_bound - T1_val;
err_low  = T1_val - T1_low_bound;

fprintf('T1 (50%% CI): %.2f (+%.2f / -%.2f) us\n', T1_val, err_high, err_low);


figure('Color', 'w'); 
hold on;

hData = plot(t, norm_contrast_cum, 'k-o', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

% Error Bars 
hError = errorbar(t, norm_contrast_cum, norm_std_dev, 'Color', [0.5 0.5 0.5], 'LineStyle', 'none'); %norm_std_err*P95

% Plot Fit Curve

hFit = plot(curve); 
set(hFit, 'Color', '#FF8800', 'LineWidth', 4);

box on; 
grid off;
set(gca, 'LineWidth', 1.2, 'FontSize', 12, 'TickDir', 'out', 'XMinorTick', 'on', 'YMinorTick', 'on');


xlabel('\tau_r (\mus)', 'FontWeight', 'bold');
ylabel('NV coherence', 'FontWeight', 'bold');

% text(max(t)*0.5, 0.9, sprintf('T_1 = %.2f \\pm %.2f \\mus', T1_val, T1_SD), ...
%      'FontSize', 11, 'FontWeight', 'bold', 'EdgeColor', 'k', 'BackgroundColor', 'w');

text_str = sprintf('T_1 = %.2f ^{+%.2f}_{-%.2f} \\mus', T1_val, err_high, err_low);

text(max(t)*0.4, 0.95, text_str, ...
     'FontSize', 11, 'FontWeight', 'bold', 'EdgeColor', 'k', ...
     'BackgroundColor', 'w', 'Interpreter', 'tex');



hold off;

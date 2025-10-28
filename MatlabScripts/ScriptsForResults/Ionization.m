%% Run test cross-section
%
% This script calculates the Bremsstrahlung doubly differential 
% cross section with screening corrections following the Olsen--Maximon--Wergeland
% additivity rule. The atomic potential is described by a Multi-Yukawa
% potential. 
%
% Written by S. Guinchard 
% < salomon.guinchard@epfl.ch > 
% Last update (10/28/2025)

%% RUN CROSS SECTION CALCULATION

% =========================================================================
% Constants 
% =========================================================================

[~,~,~,~,~,~,re,mc2,c,alpha] = pc_dke_yp; % Physics constant

% =========================================================================
% Numerical parameters
% =========================================================================

nk = 40;                % Number of photon energy points
E0 = 4540;               % Initial energy of incident electron (keV)
Z0 = 79;                 % Atomic number of gold
theta_arr = [0,1.45,3.01,6.03];          % Observation angles in degrees
params = struct;  
params.nion = [0:7:78]; % Ionization states
smodel = struct; 
    smodel.YKorder = 1;        % Yukawa model order
k_in = linspace(10, E0-5, nk);  % Photon energies in keV
vE = kinetic_energy_to_velocity(E0);
BornCond = alpha*Z0*c/vE;

% display electron velocity and validity limits
disp(['Electron velocity: ', num2str(vE/c), '*c'])
disp(['Born approximation validity check for incoming e-: ', num2str(BornCond)])
if BornCond > 0.33
    warning('Born approximation may not hold because Z too high or E0 too low.') 
end

%% Yukawa Data

YukawaModel = load_yukawa('.', Z0, params.nion, smodel.YKorder, true);

figure

% --------------------------------------------------
% First subplot: Lambdas
% --------------------------------------------------
subplot(1,2,1)   % left panel
for i = 1:smodel.YKorder
    plot(0:1:Z0-1,YukawaModel.Lambda(:,i), 'o-', 'LineWidth', 2, ...
         'DisplayName', ['$\bar{\Lambda}_{', num2str(i), '}$']);
    hold on
end
grid on
grid minor
xlim([0,Z0-1])
set(gca, 'FontSize', 25)
xlabel('$Z_{s^{\prime}}$', 'Interpreter', 'latex')
ylabel('$\bar{\lambda}$ [a.u.]', 'Interpreter', 'latex')
legend('Location', 'northwest', 'Interpreter', 'latex')
title('Inverse Screening Lengths $\bar{\lambda}_i$', 'Interpreter', 'latex')

% --------------------------------------------------
% Second subplot: A coefficients
% --------------------------------------------------
subplot(1,2,2)   % right panel
for i = 1:smodel.YKorder
    plot(YukawaModel.A(:,i), 'o-', 'LineWidth', 2, ...
         'DisplayName', ['$A_{', num2str(i), '}$']);
    hold on
end
plot(sum(YukawaModel.A(:,1:smodel.YKorder),2), 'k-', 'LineWidth', 1, ...
     'DisplayName', '$\sum A_i$');  % Sum
grid on
grid minor
xlim([0,Z0-1])
set(gca, 'FontSize', 25)
xlabel('$Z_{s^{\prime}}$', 'Interpreter', 'latex')
ylabel('$A_i$ ', 'Interpreter', 'latex')
legend('Location', 'northwest', 'Interpreter', 'latex')
title('Weights $A_i$', 'Interpreter', 'latex')

%% COMPUTE COULOMB CORRECTIONS OR ALTERNATE SCREENING MODELS

addpath(genpath('/Users/salomonguinchard/Documents/GitHub/MatTools'))

screening_models = {'Tseng-Pratt-Botto', 'full-screening', ...
'No-screening','Thomas-Fermi','Thomas-Fermi-Kirillov', ...
'Tseng-Pratt-Avdonina-Lamoureux', ...
'Tseng-Pratt-Thomas-Fermi', 'moliere', 'Multi-Yukawa'};

% now compute Coulomb corrected cross-sections
smodel.screening_model = screening_models{9}; 
smodel.Z0 = params.nion;

%RDP parameters
RelTol = 1e-6;
nphi = NaN; 
ntheta = NaN;    
nhyp = NaN; 

integral_mode = 'simps'; % Numerical integration : trapz (trapeze) or simps (Simpson) 
use_simps =  exist('integral_mode', 'var') && strcmp(integral_mode, 'simps') && exist('simps', 'file'); 
if use_simps
    if mod(nphi,2) == 0 
        nphi = nphi+1;
        warning('nphi increased to %d to satisfy Simpson''s rule requirement.', nphi);
    end
    if mod(ntheta,2) == 0
        ntheta = ntheta+1;
     warning('ntheta increased to %d to satisfy Simpson''s rule requirement.', ntheta);
    end
end

% PARAMS FOR BHE_DKE_KK
% Memory management enforcement :
%  0 -> 5-D block numerical calculations, 
%  1 -> 4-D block numerical calculations 
flag_mem = 0;                 
if exist('flag_mem', 'var') 
    params.flag_mem = flag_mem; 
end 
params.ntheta = 751; 
params.nphi = 1501;
params.nhyp = 200; 
params.integration_method = integral_mode;


% Below
% 1st order Coulomb corrections (BHE)
% Screened cross-section (BHS-BHES)
% 2nd order Coulomb corrections (RDP)
tic
for j = 1:numel(theta_arr)
   for ik = 1:nk
        % analytical BH DDCS 
        %-------------------
        [cs_BHE(j,1,ik), cs_BH(j,1,ik), ~, ~, ~, ~] = bhe_dke_yp(E0, k_in(ik), cos(deg2rad(theta_arr(j))), Z0);
        
        % RDP DDCS 
        %------------------
        [cs_RDP(j,1,ik),cs_EH(j,1,ik),cs_RDP2(j,1,ik)] = rdp_dke_jd(E0,k_in(ik),cos(deg2rad(theta_arr(j))), Z0, nhyp, nphi, ntheta, integral_mode,RelTol*cs_BHE(j,1,1),RelTol); 
   
        % DDCS with alternate screening models 
        %-------------------------------------
        [cs_BHES(:,j,ik),cs_BHS(:,j,ik),~,~,~,~] = bhe_dke_yp(E0,k_in(ik),cos(theta_arr(j)*pi/180),Z0, params, smodel);
        
   end
end
parfor j = 1:numel(theta_arr)
   for ik = 1:nk
        % DDCS with alternate screening models 
        %-------------------------------------
        [cs_BHES(:,j,ik),cs_BHS(:,j,ik),~,~,~,~] = bhe_dke_yp(E0,k_in(ik),cos(theta_arr(j)*pi/180),Z0, params, smodel);
        
   end
end
toc

%% Analytical doubly differential cross section with MY

YukawaModel = load_yukawa('.',Z0, params.nion, smodel.YKorder, false); 
A = YukawaModel.A(params.nion+1,:);                  % Extract Yukawa Weights 
a = 2./alpha./YukawaModel.Lambda(params.nion+1,:);   % Yukawa screening length  
b_bar = alpha*YukawaModel.Lambda(params.nion+1, :);  % shape (nion, YKorder)
b_bar2 = b_bar.^2;
c_bar2 = zeros(size(b_bar2));
qZs = params.nion/Z0;
for ii = 1:numel(params.nion)
    c_bar2(ii,:) =  qZs(ii) * b_bar2(ii,:);
end

cs_BHS_a = zeros([numel(theta_arr), numel(params.nion), numel(k_in)]);
cs_BHS_a(1,:,:) = bhs_my_sg(E0, k_in, deg2rad(theta_arr(1)), Z0, smodel.YKorder, params.nion, A, b_bar, c_bar2);  
if numel(theta_arr) > 1
    for ii = 2:numel(theta_arr)
        cs_BHS_a(ii,:,:) = bhs_my_sg(E0, k_in, deg2rad(theta_arr(ii)), Z0, smodel.YKorder, params.nion, A, b_bar, c_bar2);  
    end
end 

%% Total doubly differential cross section
for j = 1:numel(params.nion)
    for ii=1:numel(theta_arr)
        cs_Total(ii,:) = squeeze(cs_RDP(ii,1,:)) + squeeze(cs_BHS_a(ii,j,:)) - squeeze(cs_BH(ii,1,:));
    end
end

%% PLOT IONIZED CROSS SECTIONS

% For color 
ncolors = 128;  % Number of colors in the colormap
% Define start and end RGB colors
blue = [0, 0, 1];       % 'b'
cyan = [0, 1, 1];       % 'c'
% Interpolate between blue and cyan
custom_cmap = [linspace(blue(1), cyan(1), ncolors)', ...
               linspace(blue(2), cyan(2), ncolors)', ...
               linspace(blue(3), cyan(3), ncolors)'];
% Shift and map to theta array
inward_shift = 1;
theta_idx = round(linspace(0 + inward_shift, ncolors - 1 - inward_shift, numel(params.nion))) + 1;

% Plot
theta_ind = 4; % angle to plot

expdatas = exp_FEB_results(8);
n_arr = expdatas.n;
thetas = expdatas.theta;
plot_exp = false;

figure
cmap = custom_cmap(theta_idx, :);
rainbow = colormap(custom_cmap); 
for j = theta_ind
    for ii = 1:numel(params.nion)
        semilogy(k_in / 1000,  squeeze(cs_BHS(ii,theta_ind,:)) * 1e4 / 0.511  , '-+', ...
         'linewidth', 3, ...
         'DisplayName', [smodel.screening_model,', nion=' num2str(params.nion(ii))], ...
         'Color',  rainbow(theta_idx(ii), :));
        hold on
        semilogy(k_in / 1000,  squeeze(cs_BHS_a(theta_ind,ii,:)) * 1e4 / 0.511, '-', ...
             'linewidth', 3, ...
             'DisplayName', ['Analytical, nion=' num2str(params.nion(ii))], ...
             'Color',  rainbow(theta_idx(ii), :));
    hold on 
    end
    semilogy(k_in / 1000,  squeeze(cs_BH(theta_ind,:)) * 1e4 / 0.511, '-', ...
     'linewidth', 3, ...
     'DisplayName', 'BH', ...
     'Color',  'r');
    hold on 
end
set(gca, 'FontSize', 25);
set(gca, 'YScale', 'log');
xlim([0,max(k_in)/1000])
xlabel('$k$ [MeV]', 'Interpreter', 'latex');
ylabel('$\frac{d^2\sigma}{dk \, d\Omega_k}$ [cm$^2$/sr/MeV]', 'Interpreter', 'latex');
title(['DDCS, ', '$\theta_0=$', num2str(theta_arr(theta_ind)), ', ' ...
       '$Z=$', num2str(Z0), ' $E_c=$', num2str(E0/1000), ' MeV'], 'Interpreter', 'latex');
legend('Location', 'northeast');
grid on; grid minor;

%% PLOT IONIZED CROSS SECTIONS LIN SCALE

k_shift = 1;
normalization_fact = k_in(k_shift:end)'*1e28 /511; %/Z0^2;
expdatas = exp_FEB_results(9);
n_arr = expdatas.n;
thetas = expdatas.theta;
theta_ind = 3;
fig = figure('Color','w');
hold on
plot_exp = false;
if plot_exp
    for ii = theta_ind
        [xp, idx] = sort(expdatas.sigma{theta_ind}(:,1));  % xp in MeV
        yp = expdatas.sigma{theta_ind}(idx, 2);            % cm² / MeV / sr
        [stat_err, syst_err] = errors_cross_sec(E0/1000, xp, yp);

        % Convert: cm² / (MeV·sr) × MeV × 1e-24 → barn / sr
        scale_y = 10^(-n_arr(theta_ind)) * xp * 1e27 /511 * 0.511; % / Z0^2;

        yp = yp .* scale_y;
        stat_err = stat_err.* scale_y ;%.* scale_y;
        syst_err = syst_err.* scale_y ;%.* scale_y;

        h_fill = fill([xp; flipud(xp)], [yp + syst_err; flipud(yp - syst_err)], ...
            [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
            'HandleVisibility', 'off');
        hold on
        h = errorbar(xp, yp, stat_err, 'ko', 'MarkerFaceColor', 'w', ...
            'LineWidth', 1.2, 'MarkerSize', 12, ...
            'HandleVisibility', 'off');

        hold on;
    end
end
for ii = 1:numel(params.nion)
    plot(k_in(k_shift:end)/1000, normalization_fact.* squeeze(cs_BHS(ii,1,k_shift:end)), '+-', ...
        'Color',  rainbow(theta_idx(ii), :), 'DisplayName', ['N=', num2str(params.nion(ii))], 'linewidth', 2);
     hold on 
    plot(k_in(k_shift:end)/1000, normalization_fact.* squeeze(cs_BHS_a(1,ii,k_shift:end)), ...
        '-','Color',  rainbow(theta_idx(ii), :), 'DisplayName', ['nion=' num2str(params.nion(ii))], 'linewidth', 2);
    hold on
end
plot(k_in(k_shift:end)/1000, normalization_fact.* squeeze(cs_BH(1,1,k_shift:end)), ...
        'Color',  'r', 'DisplayName', 'Bethe-Heitler-Sauter', 'linewidth', 2);
set(gca, 'FontSize', 25);
xlim([0, max(k_in(end)/1000)])
xlabel('$k$ [MeV]', 'Interpreter', 'latex');
ylabel('$k\frac{d^2\sigma}{dk \, d\Omega_k}$ [barn/sr]', 'Interpreter', 'latex');
title(['$Z=$', num2str(Z0), ', $E_c=$', num2str(E0/1000), ' MeV, ', '$\theta_0 =$', num2str(thetas(theta_ind)), '$^\circ$'], ...
        'Interpreter', 'latex');
legend('Location', 'northeast');
grid on; grid minor;


%% Plot screening effect region

% Compute ratio of screened to unscreened BH 
ratio_map = zeros(numel(theta_arr), nk);
for j = 1:numel(theta_arr)
    ratio_map(j,:) = squeeze(cs_BH(j,1,:)) ./ squeeze(cs_BHS(j,1,:) );  % 2nd vs 1st ion state
    %ratio_map(j,:) = squeeze(cs_BHS(j,1,:)) ./ squeeze(cs_Total(j,:)' );  % 2nd vs 1st ion state
    ratio_map(ratio_map<0)=0;
end

gridres = 2000;
figure
% Mesh and interpolation
[X, Y] = meshgrid(k_in/1000, theta_arr);
k_fine = linspace(min(k_in)/1000, max(k_in)/1000, gridres);
theta_fine = linspace(min(theta_arr), max(theta_arr), gridres);
[Xq, Yq] = meshgrid(k_fine, theta_fine);
Zq = interp2(X, Y, log10(ratio_map), Xq, Yq, 'spline');

% Plot
imagesc(k_fine, theta_fine, Zq);
set(gca, 'YDir', 'normal');
colorbar;

% Labels and title with LaTeX
xlabel('$k~[\mathrm{MeV}]$', 'Interpreter', 'latex');
ylabel('$\theta~[^\circ]$', 'Interpreter', 'latex');  % Use degree symbol
title('$\log_{10} \left( \frac{\sigma_{\mathrm{screened}}}{\sigma_{\mathrm{unscreened}}} \right)$', ...
      'Interpreter', 'latex');

% Aesthetics
set(gca, 'FontSize', 25)

%
% Interpolation Grid
[X, Y] = meshgrid(k_in/1000, theta_arr);
k_fine = linspace(min(k_in)/1000, max(k_in)/1000, gridres);
theta_fine = linspace(min(theta_arr), max(theta_arr), gridres);
[Xq, Yq] = meshgrid(k_fine, theta_fine);
Zq = interp2(X, Y, log10(ratio_map), Xq, Yq, 'spline');

% Plotting
figure
tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

% === Left: Full k range ===
nexttile;
imagesc(k_fine, theta_fine, Zq);
set(gca, 'YDir', 'normal');
xlabel('$k~[\mathrm{MeV}]$', 'Interpreter', 'latex');
ylabel('$\theta~[^\circ]$', 'Interpreter', 'latex');
title('Full $k$ Range', 'Interpreter', 'latex');
colorbar;
set(gca, 'FontSize', 22);

% === Right: Zoomed-in k range ===
nexttile;
imagesc(k_fine, theta_fine, Zq);
set(gca, 'YDir', 'normal');
xlabel('$k~[\mathrm{MeV}]$', 'Interpreter', 'latex');
title('$k \in [k_{\min}, 2]~\mathrm{MeV}$', 'Interpreter', 'latex');
colorbar;
xlim([min(k_in)/1000, 0.8]);  % Zoom to k_min to 2 MeV
set(gca, 'FontSize', 22);

%% Plot Coulomb corrections region

% Compute ratio of screened to unscreened BH 
ratio_map = zeros(numel(theta_arr), nk);
for j = 1:numel(theta_arr)
    ratio_map_RDP(j,:) = squeeze(cs_RDP(j,1,:)) ./ squeeze(cs_BH(j,1,:));  % 2nd vs 1st ion state
    ratio_map_RDP(ratio_map_RDP<0)=0;
end

gridres=4000;

figure
% Mesh and interpolation
[X, Y] = meshgrid(k_in/1000, theta_arr);
k_fine = linspace(min(k_in)/1000, max(k_in)/1000, gridres);
theta_fine = linspace(min(theta_arr), max(theta_arr), gridres);
[Xq, Yq] = meshgrid(k_fine, theta_fine);
ZqRDP = interp2(X, Y, log10(ratio_map_RDP), Xq, Yq, 'spline');

% Plot
imagesc(k_fine, theta_fine, ZqRDP);
set(gca, 'YDir', 'normal');
colorbar;

% Labels and title with LaTeX
xlabel('$k~[\mathrm{MeV}]$', 'Interpreter', 'latex');
ylabel('$\theta~[^\circ]$', 'Interpreter', 'latex');  % Use degree symbol
title('$\log_{10} \left( \frac{\sigma_{\mathrm{screened}}}{\sigma_{\mathrm{unscreened}}} \right)$', ...
      'Interpreter', 'latex');

% Aesthetics
set(gca, 'FontSize', 25)

%
% Interpolation Grid
[X, Y] = meshgrid(k_in/1000, theta_arr);
k_fine = linspace(min(k_in)/1000, max(k_in)/1000, gridres);
theta_fine = linspace(min(theta_arr), max(theta_arr), gridres);
[Xq, Yq] = meshgrid(k_fine, theta_fine);
ZqRDP = interp2(X, Y, log10(ratio_map_RDP), Xq, Yq, 'spline');

% Plotting
figure
tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

% === Left: Full k range ===
nexttile;
imagesc(k_fine, theta_fine, ZqRDP);
set(gca, 'YDir', 'normal');
xlabel('$k~[\mathrm{MeV}]$', 'Interpreter', 'latex');
ylabel('$\theta~[^\circ]$', 'Interpreter', 'latex');
title('Full $k$ Range', 'Interpreter', 'latex');
colorbar;
set(gca, 'FontSize', 22);

% === Right: Zoomed-in k range ===
nexttile;
imagesc(k_fine, theta_fine, ZqRDP);
set(gca, 'YDir', 'normal');
xlabel('$k~[\mathrm{MeV}]$', 'Interpreter', 'latex');
title('$k \in [k_{max}-0.1, k_{max}]~\mathrm{MeV}$', 'Interpreter', 'latex');
colorbar;
xlim([(max(k_in)-100)/1000, (max(k_in))/1000]);  % Zoom to k_min to 2 MeV
set(gca, 'FontSize', 22);

%% Plot screening region and Coulomb region

colors = 'jet';

figure
tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

% === Left: Screening ratio ===
ax1 = nexttile;
imagesc(k_fine, theta_fine, Zq);  
set(ax1, 'YDir', 'normal');
colormap(ax1, colors);
set(gca, 'YDir', 'normal');
xlabel('$k~[\mathrm{MeV}]$', 'Interpreter', 'latex');
ylabel('$\theta_0~[^\circ]$', 'Interpreter', 'latex');
title('Screening effects $\log_{10}(\frac{\sigma_{BH}}{\sigma_{BHS}})$', 'Interpreter', 'latex');
colorbar;
xlim([min(k_in)/1000, 0.8]);
set(gca, 'FontSize', 22);

% === Right: RDP ratio ===
ax2 = nexttile;
imagesc(k_fine, theta_fine, ZqRDP);
set(ax1, 'YDir', 'normal');
colormap(ax2, colors);
set(gca, 'YDir', 'normal');
xlabel('$k~[\mathrm{MeV}]$', 'Interpreter', 'latex');
title('Coulomb effects: $\log_{10}(\frac{\sigma_{RDP}}{\sigma_{BH}})$', 'Interpreter', 'latex');
colorbar;
xlim([(max(k_in)-100)/1000, (max(k_in))/1000]);
set(gca, 'FontSize', 22);


%% IONIZATION DEPENDENCE: Cross section vs nion

fixed_k_idx = round(numel(k_in)*[1/10,1/5,1/2,4/5]);  % pick mid-energy
dsigma = zeros(numel(fixed_k_idx), numel(params.nion));
figure;
theta_ind = 1;
for j = 1:numel(fixed_k_idx) 
        normalization_fact = k_in(fixed_k_idx(j))'*1e28 /511; %/Z0^2;
        dsigma(j,:) = normalization_fact.*squeeze(cs_BHS_a(theta_ind,:,fixed_k_idx(j)));
        hold on
        semilogy(params.nion, 10^(n_arr(i))*squeeze(cs_BHS_a(theta_ind,:,fixed_k_idx(j))) * 1e4 / 0.511, 'o-', 'LineWidth', 2, 'DisplayName', ['\theta=', num2str(theta_arr(i)), '°', ', k = ' num2str(k_in(fixed_k_idx(j))/1000)], 'linewidth', 2)
end
xlabel('$Z_{s^\prime,s}$', 'Interpreter', 'latex'); 
xlim([0,max(Z0-params.nion)])
ylabel('$k\frac{d^2\sigma}{dk \, d\Omega_k}$ [barn/sr]', 'Interpreter', 'latex');
title(['Z=', num2str(Z0),  ', $d^2\sigma/dkd\Omega_k$ vs. nion at k = ' num2str(k_in(fixed_k_idx)/1000) ' MeV'], 'Interpreter', 'latex');
legend
grid on;
set(gca, 'FontSize', 25);

colors = [1 0 0 
          0 0 1
          0.3010 0.7450 0.9330 
          0 0 0];
plot_four_curves(params.nion, dsigma(1,:), dsigma(2,:), dsigma(3,:), dsigma(4,:), CMap=colors, Karr=k_in(fixed_k_idx))

% This script calculates the Bremsstrahlung doubly differential 
% cross section with screening corrections following the Olsen--Maximon--Wergeland
% additivity rule. The atomic potential is described by a Multi-Yukawa
% potential. The total model is then compared to experiments.
%
% Written by S.Guinchard and J.Decker 
% < salomon.guinchard@epfl.ch > < joan.decker@epfl.ch >
% Last update (10/28/2025)

clear
close all
% =========================================================================
% Constants 
% =========================================================================
alpha = 1/137;   % fine-structure constant
re = 2.817e-15;  % m, classical electron radius
keV = 1.602e-16; % J, keV conversion
me = 9.1e-31;    % kg, electron mass
mbarn = 1e-31;   % m^2, 1 millibarn
c = 3e8;         % m/s, speed of light
a0 = 5.29e-11;   % m, Bohr radius

% =========================================================================
% Numerical parameters
% =========================================================================

nk = 200;                % Number of photon energy points

ind_exp = [20,9,17,23];     % Exp data index
theta_sel = {...% selection of theta0 for calculation and display (for each experimental dataset). 
[],...   %  Leave empty to select all  
[3,7],...% 10° and 60°
[3,7],...% 10° and 60°
[3,7]};  % 10° and 60°

for ind=1:numel(ind_exp)

    expdatas = exp_FEB_results(ind_exp(ind));

    if isempty(theta_sel{ind})
        tmask = 1:length(expdatas.theta);
    else
        tmask = theta_sel{ind};
    end

    E0 = expdatas.E0*1000;   % Initial energy of incident electron (keV)
    Z0 = expdatas.Znucleus;  % Atomic number
    
    n_arr = expdatas.n(tmask);
    theta_arr = expdatas.theta(tmask);
    sigma = expdatas.sigma(tmask);
    plot_exp = true;
    nj = numel(theta_arr);
    
    params = struct;  
    params.nion = [0];               % Ionization states
    smodel = struct; 
        smodel.YKorder = 3;          % Yukawa model order
    k_in = linspace(10, E0 - 5, nk);  % Photon energies in keV
    vE = kinetic_energy_to_velocity(E0);
    BornCond = alpha*Z0*c/vE;
    
    % display electron velocity and validity limits
    disp(['Electron velocity: ', num2str(vE/c), '*c'])
    disp(['Born approximation validity check for incoming e-: ', num2str(BornCond)])
    if BornCond > 0.33
        warning('Born approximation may not hold because Z too high or E0 too low.') 
    end

    screening_models = {'Tseng-Pratt-Botto', 'full-screening', ...
    'No-screening','Thomas-Fermi','Thomas-Fermi-Kirillov', ...
    'Tseng-Pratt-Avdonina-Lamoureux', ...
    'Tseng-Pratt-Thomas-Fermi'  };
    
    % now compute Coulomb corrected cross-sections
    smodel.screening_model = screening_models{1};  % Tseng-Pratt-Botto
    
    nphi = NaN;%751; 
    ntheta = NaN;%1501;    
    nhyp = NaN;%200; 

    RelTol = 1e-6;

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
    flag_mem = 1;                 % Memory management enforcement : 0 -> 5-D block numerical calculations, 1 -> 4-D block numerical calculations 
    if exist('flag_mem', 'var') 
        params.flag_mem = flag_mem; 
    end 
    params.ntheta = ntheta; 
    params.integral = integral_mode;
    
    % Below
    % 1st order Coulomb corrections (BHE)
    % Screened cross-section (BHS-BHES)
    % 2nd order Coulomb corrections (RDP)
    
    cs_RDP = zeros(nj,1,nk);
    cs_EH = zeros(nj,1,nk);
    cs_RDP2 = zeros(nj,1,nk);

    for j = 1:nj
        for ik = 1:nk
            [cs_BHE(j,1,ik), cs_BH(j,1,ik), ~, ~, ~, ~] = bhe_dke_yp(E0, k_in(ik), cos(deg2rad(theta_arr(j))), Z0);
            [cs_RDP(j,1,ik),cs_EH(j,1,ik),cs_RDP2(j,1,ik)] = rdp_dke_jd(E0,k_in(ik),cos(deg2rad(theta_arr(j))), Z0, nhyp, nphi, ntheta, integral_mode, RelTol*cs_BHE(j,1,1),RelTol); 
          % [se_BHES(:,1,ik),se_BHS(:,1,ik),~,~,~,~] = bhe_dke_yp(E0,k_in(ik),cos(theta_arr*pi/180),Z0, params, smodel);
            disp(['processed case j=',num2str(j),'/',num2str(nj),'; k=',num2str(ik),'/',num2str(nk)])            
        end
    end

    % Below
    % Screened Bethe-Heitler (Born) cross-section
    % Model: Multi-Yukawa
    YukawaModel = load_yukawa('.',Z0, params.nion, smodel.YKorder, false); 
    A = YukawaModel.A(params.nion+1,:);                  % Extract Yukawa Weights 
    a = 2./alpha./YukawaModel.Lambda(params.nion+1,:);   % Yukawa screening length  
    b_bar = alpha*YukawaModel.Lambda(params.nion+1, :); % shape (nion, YKorder)
    b_bar2 = b_bar.^2;
    c_bar2 = zeros(size(b_bar2));
    qZs = params.nion/Z0;
    for ii = 1:numel(params.nion)
        c_bar2(ii,:) =  qZs(ii) * b_bar2(ii,:);
    end
    
    cs_BHS_a = zeros([nj, numel(params.nion), nk]);
    cs_BHS_a(1,:,:) = bhs_my_sg(E0, k_in, deg2rad(theta_arr(1)), Z0, smodel.YKorder, params.nion, A, b_bar, c_bar2);  
    if nj > 1
        for ii = 2:nj
            cs_BHS_a(ii,:,:) = bhs_my_sg(E0, k_in, deg2rad(theta_arr(ii)), Z0, smodel.YKorder, params.nion, A, b_bar, c_bar2);  
        end
    end 
%%
    % total cross-section
    
    for ii=1:nj
        cs_Total(ii,:) = squeeze(cs_RDP(ii,1,:)) + squeeze(cs_BHS_a(ii,1,:))-squeeze(cs_BH(ii,1,:));
    end
    
    % COMPONENT ANALYSIS (Individual)

    path = 'Fig_test';

    for theta_ind = 1:nj
        k_MeV = k_in / 1000;
        k_shift=1; % shift to avoid infrared divergence plot
        % color parameters for plotting
        
        ncolors = 128;
        inward_shift = 30;
        rainbow = spring(ncolors); 
        theta_norm = (theta_arr - min(theta_arr)) / (max(theta_arr) - min(theta_arr));  
        theta_idx = round(linspace(0+inward_shift, ncolors-1-inward_shift, nj)) + 1;
        cmap = spring(nj);
        
            
            fig = figure;
            semilogy(k_in/1000, 10^(n_arr(theta_ind)).*squeeze(cs_RDP(theta_ind,1,:)) * 1e4 / 0.511, 'b', 'DisplayName', 'RDP', 'linewidth', 2);
            hold on;
            semilogy(k_in/1000, 10^(n_arr(theta_ind)).*squeeze(cs_BHS_a(theta_ind,1,:)) * 1e4 / 0.511, 'g', 'DisplayName', 'Screened BH', 'linewidth', 2);
            semilogy(k_in/1000, 10^(n_arr(theta_ind)).*squeeze(cs_BH(theta_ind,1,:))  * 1e4 / 0.511, 'r', 'DisplayName', 'Unscreened BH', 'linewidth', 2);
            semilogy(k_in/1000, 10^(n_arr(theta_ind)).*squeeze(real(cs_EH(theta_ind,:))) * 1e4 / 0.511, 'c', 'DisplayName', 'EH', 'linewidth', 2);
            %
            semilogy(k_in/1000, 10^(n_arr(theta_ind)).*squeeze(cs_Total(theta_ind,:)) * 1e4 / 0.511, 'k--', 'DisplayName', 'Total', 'linewidth', 2);

                [xp, idx] = sort(sigma{theta_ind}(:,1)); 
                yp = sigma{theta_ind}(idx, 2);
                [stat_err, syst_err] = errors_cross_sec(E0/1000, xp, yp);
                h_fill =  fill([xp; flipud(xp)], [yp + syst_err; flipud(yp - syst_err)], ...
                    [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
                    'HandleVisibility', 'off');
                hold on 
                h = errorbar(xp, yp, stat_err, 'ko', 'MarkerFaceColor', 'w', ...
                    'LineWidth', 1.2, 'MarkerSize', 12, ...
                    'HandleVisibility', 'off');

            set(gca, 'FontSize', 30);
            set(gca,'TickLabelInterpreter','latex')
            set(gca, 'YScale', 'log');
            grid on; grid minor;
            xlim([0, max(k_MeV)])
            xlabel('$k$ [MeV]', 'Interpreter', 'latex');
            ylabel('$10^n\frac{d^2\sigma}{dk \, d\Omega_k}$ [cm$^2$/sr/MeV]', 'Interpreter', 'latex');
            title([['$\theta_0 =$', num2str(theta_arr(theta_ind))] ...
                   ', $Z=$', num2str(Z0), ' $E_c=$', num2str(E0/1000), ' MeV'], 'Interpreter', 'latex');
            legend('Location', 'northeast', 'interpreter', 'latex');
            savefig(fig,[path,'/dsigma_all_Z=' num2str(Z0),'E=',num2str(E0/1000), 'MeV_theta',num2str(theta_arr(theta_ind)),'_log.fig']);
            close(fig);
    
    
            normalization_fact = k_in(k_shift:end)'*1e28 /511; %/Z0^2;
            path = 'Fig_test';
            fig = figure;
            hold on
                [xp, idx] = sort(sigma{theta_ind}(:,1));  % xp in MeV
                yp = sigma{theta_ind}(idx, 2);            % cm² / MeV / sr
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


             plot(k_in(k_shift:end)/1000, normalization_fact.* squeeze(cs_BH(theta_ind,1,k_shift:end)), 'r-', 'DisplayName', 'BH', 'linewidth', 2);
             hold on;
             plot(k_in(k_shift:end)/1000, normalization_fact.* squeeze(cs_BHS_a(theta_ind,1,k_shift:end)),'k--', 'DisplayName', 'BH$^{s}$', 'linewidth', 1);
             plot(k_in(k_shift:end)/1000, normalization_fact.* squeeze(cs_RDP(theta_ind,1,k_shift:end)), 'b-', 'DisplayName', 'RDP', 'linewidth', 1);
             plot(k_in(k_shift:end)/1000, normalization_fact'.* squeeze(cs_Total(theta_ind,k_shift:end)), '--','Color',[0,0.5,0], 'DisplayName', 'Total', 'linewidth', 2);
        
            set(gca, 'FontSize', 30);
            set(gca,'TickLabelInterpreter','latex')
            xlim([0, max(k_MeV)])
            xlabel('$k$ [MeV]', 'Interpreter', 'latex');
            ylabel('$k\frac{d^2\sigma}{dk \, d\Omega_k}$ [barn/sr]', 'Interpreter', 'latex');
            title(['$Z=$', num2str(Z0), ', $E_c=$', num2str(E0/1000), ' MeV, ', '$\theta_0 =$', num2str(theta_arr(theta_ind)), '$^\circ$'], ...
                    'Interpreter', 'latex');
            legend('Location', 'northeast', 'Interpreter', 'latex');
            grid on; grid minor;
            figname = ['dsigma_all_Z=' num2str(Z0),'E=',num2str(E0/1000), 'MeV_theta',num2str(theta_arr(theta_ind)),'_lin'];
            savefig(fig,[path,filesep,figname,'.fig']);
            exportgraphics(fig,[path,filesep,figname,'.eps'])
            % close(fig);
    end

end


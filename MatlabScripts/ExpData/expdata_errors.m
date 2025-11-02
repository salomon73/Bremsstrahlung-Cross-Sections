%% expdata_erros
%
% This scripts enables to plot the experimental data from various sources
% for bremsstrahlung differential cross sections with systematic and
% statistical errors determined from the source method. 
% 
%       
%   Dependencies: 
%
%       - errors_cross_sec.m    %  Sources of errors from exp papers
%       - exp_FEB_results.m     %  Experimental data
%
%
%   By S. Guinchard <salomon.guinchard@epfl.ch> 
%   Last update (28/10/2025)


expdatas = exp_FEB_results(19);
n_arr = expdatas.n;
theta_arr = expdatas.theta;
figure;
for ii =1:numel(n_arr)
    [xp, idx] = sort(expdatas.sigma{ii}(:,1)); 
     yp = expdatas.sigma{ii}(idx, 2);
    
    % Normalized photon energy
    E0 = expdatas.E0; % MeV
    k_over_E0 = xp / E0;
    
    [stat_err, syst_err] = errors_cross_sec(E0,xp,yp);
  
    hold on

    y_upper = yp + stat_err + syst_err;
    y_lower = yp - stat_err - syst_err;
    x_contour = [xp; flipud(xp)];           % goes forward then backward
    y_contour = [y_upper; flipud(y_lower)]; % top then bottom
    
    h_fill = patch(x_contour, y_contour, [0.8 0.8 0.8], ...
               'EdgeColor', 'k', ...    
                'linewidth', 1, ...
               'FaceAlpha', 0, ...
               'DisplayName', ['Systematic, \theta=' num2str(theta_arr(ii)) ', n=' num2str(n_arr(ii))] ...
               ,'HandleVisibility', 'off');

    % Plot experimental points with statistical error bars
    h_err = errorbar(xp, yp, stat_err, 'o', 'MarkerFaceColor', 'w', ...
    'LineWidth', 1.2, 'MarkerSize', 8, ...
    'DisplayName', ['statistical, \theta=' num2str(theta_arr(ii)) ', n=' num2str(n_arr(ii))]);

    xlim([0,E0]);
    set(gca, 'YScale', 'log');
    set(gca, 'FontSize', 25);
    xlabel('$k$ [MeV]', 'Interpreter', 'latex');
    ylabel('$10^n\frac{d^2\sigma}{dk \, d\Omega_k}$ [cm$^2$/sr/MeV]', 'Interpreter', 'latex');
    title([num2str(expdatas.elements),' $E_0$=', num2str(E0), ' MeV'], 'Interpreter', 'latex');
    grid on;
    box on;
    legend('Location', 'northeast');

end


%% Linear-scale plot

% choose theta0 angle to plot 
theta_ind =1;

figure;
hold on

    % Extract experimental data
    [xp, idx] = sort(expdatas.sigma{theta_ind}(:,1));
    yp = expdatas.sigma{theta_ind}(idx,2);

    xp = xp(:);
    yp = yp(:);

    % Errors
    [stat_err, syst_err] = errors_cross_sec(E0, xp, yp);
    stat_err = stat_err(:);
    syst_err = syst_err(:);

    % Scaling including photon energy
    scale_y = 10^(-n_arr(theta_ind)) * xp * 1e27 / 511 * 0.511;
    yp_scaled = yp .* scale_y;
    stat_err_scaled = stat_err .* scale_y;
    syst_err_scaled = syst_err .* scale_y;

    % Patch for total error
    y_upper = yp_scaled + stat_err_scaled + syst_err_scaled;
    y_lower = yp_scaled - stat_err_scaled - syst_err_scaled;
    x_contour = [xp; flipud(xp)];        % goes forward then backward
    y_contour = [y_upper; flipud(y_lower)]; % top then bottom

    h_fill = patch(x_contour, y_contour, [0.8 0.8 0.8], ...
               'EdgeColor', 'k', ...    % contour line color
                'linewidth', 1, ...
               'FaceAlpha', 0, ...
               'EdgeAlpha', 0.7, ...
               'DisplayName', ['Systematic, \theta=' num2str(theta_arr(theta_ind)) ', n=' num2str(n_arr(theta_ind))] ...
               ,'HandleVisibility', 'off');

    hold on;

    % Statistical error bars
    errorbar(xp, yp_scaled, stat_err_scaled, 'ko', 'MarkerFaceColor','w', ...
        'LineWidth',1.2, 'MarkerSize',8, 'HandleVisibility','on', ...
        'DisplayName',['$\theta_{0}=$', num2str(theta_arr(theta_ind)),'$^{\circ}$']);

xlim([0, 1.1*max(xp)]);
xlabel('$k$ [MeV]', 'Interpreter','latex');
ylabel('$k\frac{d^2\sigma}{dk \, d\Omega_k}$ [barn/sr]', 'Interpreter','latex');
title(['$E_0=$',num2str(E0),' MeV'],'Interpreter','latex');
legend('Location','northeast', 'interpreter', 'latex');
grid on; box on; set(gca,'FontSize',25);

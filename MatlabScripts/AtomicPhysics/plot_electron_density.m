%  This scripts enables the plotting of electron densities
%  obtained from the DHFS-DFT code Gaussian for the Multi-Yukawa
%  as well as alternate atomic models. 
%  
%  Refs: 
%       - Screened Thin-Target Bremsstrahlung with Partially-Ionized High-Z
%           Species, Guinchard, Peysson, Decker, 2025.
%       - Approximate atomic models for fast computation of the Fokker–Planck 
%           equation in fusion plasmas with high-Z impurities and suprathermal electrons
%           Walkowiak, Jardin et al.
%
%   By S. Guinchard <salomon.guinchard@epfl.ch> 
%   Last update (28/10/2025)
%% Plot electronic densities
load('dft_79_Au_gold.mat');

figure
ax = gca;
i_ref_ind =  1+[0:14:78];
i_all = 1:79;
j = setdiff(i_all, i_ref_ind);
for i=j
    h = semilogy(dft.r,dft.ne(:,i), '-', 'linewidth', 2,'HandleVisibility','off'); 
    h.Color = [0.8 0.2 0.1 0.2];
    hold on
end
nLines = numel(i_ref_ind);
fullMap = highContrastBWR(256);
nRows = size(fullMap,1);
idx = round(linspace(1, nRows, nLines));
colors = fullMap(idx,:);
set(ax,'ColorOrder',colors,'NextPlot','add')

for k = 1:nLines
    i=i_ref_ind(k);
    semilogy(dft.r,dft.ne(:,i), '-', 'linewidth', 3.5, 'Color', colors(k,:) ...
        , 'displayname',['$Z_{s^{\prime}}=$',num2str(i-1)]) 
    hold on
end

ylim([5e-4, 5e4])
yticks([5e-4, 1e-2, 1e-0, 1e2, 1e4])
yticklabels({'$5\cdot 10^{-4}$', '$10^{-2}$', '$10^{0}$', '$10^{2}$', '$10^{4}$'})
grid on 
grid minor
set(gca, 'fontsize', 30)
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$r_{s^{\prime}}$ [$a_0$]','interpreter', 'latex')
ylabel('$n_{s^{\prime}}$ [$a_0^{-3}$]','interpreter', 'latex')
legend('location', 'northeast', 'interpreter', 'latex')


%% Radius of intersection of electron densities
figure
ax = gca;
hold on
grid on
grid minor

nZ = size(dft.ne, 2);        % number of charge states or configurations
re_max = nan(nZ);            % max intersection radii
re_min = nan(nZ);            % min intersection radii (optional)

% Threshold below which density is considered numerically negligible
ne_thresh = 1e-5;
i_ref_ind =  1+[0:14:78];

nLines = numel(i_ref_ind);
fullMap = highContrastBWR(256);
nRows = size(fullMap,1);
idx = round(linspace(1, nRows, nLines));
colors = fullMap(idx,:);
set(ax,'ColorOrder',colors,'NextPlot','add')

for iref = i_ref_ind
    % Identify a reasonable cutoff radius where density ~ 0
    idx_cut = find(dft.ne(:, iref) < ne_thresh, 1, 'first');
    if isempty(idx_cut)
        r_cut = max(dft.r);
    else
        r_cut = dft.r(idx_cut);
    end

    for ii = iref:nZ
        if ii == iref
            continue
        end

        % Mask out the tail region (below threshold)
        mask = (dft.r <= r_cut) & ...
               (dft.ne(:, iref) > ne_thresh) & ...
               (dft.ne(:, ii) > ne_thresh);

        if nnz(mask) < 2
            continue
        end

        % Find intersections between reference and current curve
        [r0, ~] = intersections(dft.r(mask), dft.ne(mask, iref), ...
                                dft.r(mask), dft.ne(mask, ii));

        % Keep only physically meaningful intersections (below cutoff)
        r0_valid = r0(r0 <= r_cut);

        if ~isempty(r0_valid)
            re_min(iref, ii) = min(r0_valid);
            re_max(iref, ii) = max(r0_valid);
        end

    end
    % Plot the physically relevant (outermost) intersection
    plot(dft.Z0, re_max(iref, :), 'o','MarkerEdgeColor', [0.1,0.1,0.1], ...
       'MarkerFaceColor', colors(find(i_ref_ind==iref),:), ...
       'MarkerSize', 8, 'linewidth', 1, 'displayname',['$Z_{s^{\prime}}=$',num2str(iref-1)])

end
xlim([0,78])
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 25)
xlabel('$Z_{s^{\prime}}$', 'Interpreter', 'latex', 'FontSize', 27.5)
ylabel('$\max(r_{\mathrm{cross}})$ [$a_0$]', ...
       'Interpreter', 'latex', 'FontSize', 27.5)
title('Outermost intersection radius vs. ionization state', ...
     'Interpreter', 'latex')
legend('location', 'northeast', 'interpreter', 'latex')

%% Alternate models 

% Thomas Fermi: kk = 1
% Thomas Fermi Kirillov: kk = 2
% Avdonina Lamoureux: kk = 3
% Tseng Pratt Botto: kk = 4
% Tseng Pratt Thomas Fermi: kk=5

figure
ax = gca;
i_ref_ind =  1+[0:14:78];
i_all = 1:79;
j = setdiff(i_all, i_ref_ind);
kk = 1;
for i=j
    h = semilogy(dft.ionradius{1,kk}.rn{i},dft.ionradius{1,kk}.dens{i}, '-', 'linewidth', 2,'HandleVisibility','off'); 
    h.Color = [0.8 0.2 0.1 0.2];
    hold on
end
nLines = numel(i_ref_ind);
fullMap = highContrastBWR(256);
nRows = size(fullMap,1);
idx = round(linspace(1, nRows, nLines));
colors = fullMap(idx,:);
set(ax,'ColorOrder',colors,'NextPlot','add')

for k = 1:nLines
    i=i_ref_ind(k);
    semilogy(dft.ionradius{1,kk}.rn{i},dft.ionradius{1,kk}.dens{i}, '-', 'linewidth', 3.5, 'Color', colors(k,:) ...
        , 'displayname',['$Z_{s^{\prime}}=$',num2str(i-1)]) 
    hold on
end

ylim([5e-4, 5e4])
yticks([5e-4, 1e-2, 1e-0, 1e2, 1e4])
yticklabels({'$5\cdot 10^{-4}$', '$10^{-2}$', '$10^{0}$', '$10^{2}$', '$10^{4}$'})
grid on 
grid minor
set(gca, 'fontsize', 30)
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$r_{s^{\prime}}$ [$a_0$]','interpreter', 'latex')
ylabel('$n_{s^{\prime}}$ [$a_0^{-3}$]','interpreter', 'latex')
legend('location', 'northeast', 'interpreter', 'latex')

figure
ax = gca;
hold on
grid on
grid minor

nZ = size(dft.ne, 2);        % number of charge states or configurations
re_max = nan(nZ);            % max intersection radii
re_min = nan(nZ);            % min intersection radii (optional)

% Threshold below which density is considered numerically negligible
ne_thresh = 1e-5;
i_ref_ind =  1+[0:14:78];

nLines = numel(i_ref_ind);
fullMap = highContrastBWR(256);
nRows = size(fullMap,1);
idx = round(linspace(1, nRows, nLines));
colors = fullMap(idx,:);
set(ax,'ColorOrder',colors,'NextPlot','add')

for iref = i_ref_ind
    % Identify a reasonable cutoff radius where density ~ 0
    idx_cut = find(dft.ionradius{1,kk}.dens{1,iref} < ne_thresh, 1, 'first');
    if isempty(idx_cut)
        r_cut = max(dft.ionradius{1,kk}.rn{1,iref});
    else
        r_cut = dft.ionradius{1,kk}.rn{1,iref}(idx_cut);
    end

    for ii = 1:nZ
        if ii == iref
            continue
        end

        % Determine the common mask for the two curves
        rn_ref = dft.ionradius{1,kk}.rn{1,iref};
        rn_ii  = dft.ionradius{1,kk}.rn{1,ii};
        dens_ref = dft.ionradius{1,kk}.dens{1,iref};
        dens_ii  = dft.ionradius{1,kk}.dens{1,ii};

        % Interpolate the smaller grid onto the larger one
        if numel(rn_ref) < numel(rn_ii)
            rn_common = rn_ii;
            dens_ref = interp1(rn_ref, dens_ref, rn_common, 'pchip');
            dens_ii  = dens_ii;
        else
            rn_common = rn_ref;
            dens_ii  = interp1(rn_ii, dens_ii, rn_common, 'pchip');
            dens_ref = dens_ref;
        end

        % Mask out the tail region (below threshold)
        mask = (rn_common <= r_cut) & (dens_ref > ne_thresh) & (dens_ii > ne_thresh);
        if nnz(mask) < 2
            continue
        end

        % Find intersections between reference and current curve
        [r0, ~] = intersections(rn_common(mask), dens_ref(mask), rn_common(mask), dens_ii(mask));

        % Keep only physically meaningful intersections (below cutoff)
        r0_valid = r0(r0 <= r_cut);

        if ~isempty(r0_valid)
            re_min(iref, ii) = min(r0_valid);
            re_max(iref, ii) = max(r0_valid);
        end
    end

    % Plot the physically relevant (outermost) intersection
    plot(dft.Z0, re_max(iref, :), 'o','MarkerEdgeColor', [0.1,0.1,0.1], ...
       'MarkerFaceColor', colors(find(i_ref_ind==iref),:), ...
       'MarkerSize', 8, 'linewidth', 1, 'displayname',['$Z_{s^{\prime}}=$',num2str(iref-1)])
end

xlim([0,78])
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 25)
xlabel('$Z_{s^{\prime}}$', 'Interpreter', 'latex', 'FontSize', 27.5)
ylabel('$\max(r_{\mathrm{cross}})$ [$a_0$]', ...
       'Interpreter', 'latex', 'FontSize', 27.5)
%title('Outermost intersection radius vs. ionization state', ...
%      'Interpreter', 'latex')
legend('location', 'northeast', 'interpreter', 'latex')

%% Colormaps

function cmap = redBiasedBWR(n)
    if nargin < 1, n = 256; end

    % Define control colors (RGB)
    deepBlue  = [0.230, 0.299, 0.754];
    lightBlue = [0.600, 0.760, 0.930];
    middleBlueBeige = [0.85, 0.7, 0.73];
    middleRedBeige  = [0.85, 0.7, 0.73];
    lightRed  = [0.980, 0.650, 0.470];
    deepRed   = [0.700, 0.050, 0.050];

    % Split number of colors
    nBlue  = round(0.4*n);
    nMiddle= round(0.1*n);
    nRed   = n - nBlue - nMiddle;

    % Interpolate colors
    bluePart   = interp1([1 2 3],[deepBlue; lightBlue; middleBlueBeige], ...
                         linspace(1,3,nBlue));
    redPart    = interp1([1 2 3],[middleRedBeige; lightRed; deepRed], ...
                         linspace(1,3,nRed));

    % Ensure last color is exactly deepRed
    redPart(end,:) = deepRed;
    cmap = [bluePart; redPart];
end

function cmap = highContrastBWR(n)
    if nargin < 1, n = 256; end
    deepBlue = [0.1 0.25 0.8];
    midBlue  = [0.3 0.6 0.9];
    lightBlue = [0.7 0.85 0.98];
    lightRed  = [0.98 0.75 0.55];
    midRed    = [0.9 0.4 0.2];
    deepRed   = [0.6 0.0 0.05];

    nBlue = round(0.45*n);
    nRed  = n - nBlue;
    bluePart = interp1([1 2 3],[deepBlue; midBlue; lightBlue], linspace(1,3,nBlue));
    redPart  = interp1([1 2 3],[lightRed; midRed; deepRed], linspace(1,3,nRed));
    cmap = [bluePart; redPart];
end

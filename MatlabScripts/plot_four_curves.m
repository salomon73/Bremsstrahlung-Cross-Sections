function plot_four_curves(nion, sigma1, sigma2, sigma3, sigma4, varargin)
% PLOT_FOUR_CURVES  Plot 4 curves with inner & outer y-axes and custom colormap.
%
% Usage:
%   plot_four_curves(x, y1, y2, y3, y4)                   % default options
%   plot_four_curves(..., 'Mode','gradient')              % color-each-line-with-gradient
%   plot_four_curves(..., 'CMap', myCmap, 'NColors',128)  % supply custom colormap
%   plot_four_curves(..., 'FontSize', 12, 'LeftShift',0.07, 'RightShift', 0.07)
%
% Inputs:
%   nion    - x vector
%   sigma1  - left-inner y data (vector)
%   sigma2  - right-inner y data
%   sigma3  - left-outer y data (true scale)
%   sigma4  - right-outer y data (true scale)
%
% Name-value pairs:
%   'Mode'      : 'perCurve' (default) | 'gradient'
%   'CMap'      : Mx3 custom colormap (optional). If empty -> blue->cyan generated
%   'NColors'   : number of rows when generating the default colormap (default 128)
%   'FontSize'  : base font size (default 11)
%   'LeftShift' : horizontal shift for left outer ticks (default 0.06)
%   'RightShift': horizontal shift for right outer ticks (default 0.06)

% -------------------- parse inputs --------------------
p = inputParser;
addParameter(p,'Mode','perCurve', @(s) any(validatestring(s,{'perCurve','gradient'})));
addParameter(p,'CMap',[], @(x) isempty(x) || (isnumeric(x) && size(x,2)==3));
addParameter(p,'Karr',[], @(x) isempty(x) || (isnumeric(x) && size(x,2)==4));
addParameter(p,'NColors',128, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p,'FontSize',11, @(x) isnumeric(x) && isscalar(x));
addParameter(p,'LeftShift',0.06, @(x) isnumeric(x) && isscalar(x));
addParameter(p,'RightShift',0.06, @(x) isnumeric(x) && isscalar(x));

parse(p,varargin{:});
opts = p.Results;

mode = opts.Mode;
custom_cmap = opts.CMap;
ncolors = opts.NColors;
fsize = opts.FontSize;
leftShift = opts.LeftShift;
rightShift = opts.RightShift;

% basic sanity
x = nion(:)';
yLinner = sigma1(:)';
yRinner = sigma2(:)';
yLouter = sigma3(:)';
yRouter = sigma4(:)';

if ~(numel(x)==numel(yLinner) && numel(x)==numel(yRinner) && numel(x)==numel(yLouter) && numel(x)==numel(yRouter))
    error('All input vectors must have the same length.');
end

% -------------------- build custom cmap if not provided --------------------
if isempty(custom_cmap)
    blue = [0 0 1];
    cyan = [0 1 1];
    custom_cmap = [linspace(blue(1), cyan(1), ncolors)', ...
                   linspace(blue(2), cyan(2), ncolors)', ...
                   linspace(blue(3), cyan(3), ncolors)'];
else
    % if user provided 
    ncolors = size(custom_cmap,1);
end

% -------------------- style defaults (will be overwritten if mode==perCurve) ----
% default single-colors (kept only if mode == 'perCurve')
idxs = round(linspace(1, ncolors, 4));
cLinner = custom_cmap(idxs(1), :);
cRinner = custom_cmap(idxs(2), :);
cLouter = custom_cmap(idxs(3), :);
cRouter = custom_cmap(idxs(4), :);


% -------------------- prepare figure & main axes --------------------
fig = figure('Units','normalized', 'Color','w');
ax = axes(fig,'Position',[0.16 0.16 0.68 0.72], 'ActivePositionProperty','position');
hold(ax,'on');
set(ax, 'FontSize', fsize);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 30)

% -------------------- helper: linear mapping (source -> destination) ----------
mapLinear = @(y, src, dst) ( dst(1) + ( (y - src(1)) ./ max(eps,(src(2)-src(1))) ) .* (dst(2)-dst(1)) );

% -------------------- true outer limits (pad slightly) -----------------------
leftOuterTrueLim  = [min(yLouter), max(yLouter)];
padL = 0.05 * (leftOuterTrueLim(2)-leftOuterTrueLim(1));
leftOuterTrueLim = leftOuterTrueLim + [-padL padL];

rightOuterTrueLim = [min(yRouter), max(yRouter)];
padR = 0.05 * (rightOuterTrueLim(2)-rightOuterTrueLim(1));
rightOuterTrueLim = rightOuterTrueLim + [-padR padR];

% -------------------- plotting depending on mode ---------------------------
% We'll store handles in h1..h4. For 'gradient' mode we use surface so colors
% are interpolated along x. For 'perCurve' mode we use simple plot lines.

% set colormap on the axis so 'gradient' surfaces use it
set(ax,'Colormap', custom_cmap);

switch mode
    case 'perCurve'
        % left inner
        yyaxis(ax,'left');
        h1 = plot(ax, x, yLinner, ...
            'o-','MarkerEdgeColor', [0.1,0.1,0.1], ...
            'LineWidth', 0.5, 'MarkerFaceColor', cLinner,'MarkerSize', 8, ...
            'DisplayName','Left inner', 'Color', cLinner);
        ax.YAxis(1).Color = cLinner;
        ylabel(ax,'','Color',cLinner);
        
        
        % right inner
        yyaxis(ax,'right');
        h2 = plot(ax, x, yRinner, ...
            'o-','MarkerEdgeColor', [0.1,0.1,0.1], ...
            'LineWidth', 0.5, 'MarkerFaceColor', cRinner,'MarkerSize', 8, ...
            'DisplayName','Left inner', 'Color', cRinner);
        ax.YAxis(2).Color = cRinner;
        ylabel(ax,'','Color',cRinner);
       
        
        % map outer curves into main axes coordinates and plot (single color)
        yyaxis(ax,'left'); leftInnerLim = ylim(ax);
        yLouter_mapped = mapLinear(yLouter, leftOuterTrueLim, leftInnerLim);
        h3 = plot(ax, x, yLouter, ...
            'o-','MarkerEdgeColor', [0.1,0.1,0.1], ...
            'LineWidth', 0.5, 'MarkerFaceColor', cLouter,'MarkerSize', 8, ...
            'DisplayName','Left outer', 'Color', cLouter);
        
        yyaxis(ax,'right'); rightInnerLim = ylim(ax);
        yRouter_mapped = mapLinear(yRouter, rightOuterTrueLim, rightInnerLim);
        h4 = plot(ax, x, yRouter, ...
            'o-','MarkerEdgeColor', [0.1,0.1,0.1], ...
            'LineWidth', 0.5, 'MarkerFaceColor', cRouter,'MarkerSize', 8, ...
            'DisplayName','Right outer', 'Color', cRouter);
        
    case 'gradient'
        % For gradient mode we draw each curve as a 2xN surface with 'EdgeColor','interp'
        % and CData running from 1..ncolors (mapped to custom_cmap).
        cvals = linspace(1, ncolors, numel(x));
        % left inner (native)
        yyaxis(ax,'left');
        h1 = surface(ax, [x;x], [yLinner;yLinner], zeros(2,numel(x)), ...
                     repmat(cvals,2,1), ...
                     'FaceColor','none','EdgeColor','interp','LineWidth',1.6, 'DisplayName','Left inner');
        set(ax,'CLim',[1 ncolors]);
        ax.YAxis(1).Color = [0 0 0]; % keep axis color neutral or override below
        ylabel(ax,'');
        
        % right inner (native)
        yyaxis(ax,'right');
        h2 = surface(ax, [x;x], [yRinner;yRinner], zeros(2,numel(x)), ...
                     repmat(cvals,2,1), ...
                     'FaceColor','none','EdgeColor','interp','LineWidth',1.6, 'DisplayName','Right inner');
        set(ax,'CLim',[1 ncolors]);
        ax.YAxis(2).Color = [0 0.5 0];
        ylabel(ax,'');
        
        % left outer: map to left-inner coords first then plot as gradient
        yyaxis(ax,'left'); leftInnerLim = ylim(ax);
        yLouter_mapped = mapLinear(yLouter, leftOuterTrueLim, leftInnerLim);
        h3 = surface(ax, [x;x], [yLouter_mapped;yLouter_mapped], zeros(2,numel(x)), ...
                     repmat(cvals,2,1), ...
                     'FaceColor','none','EdgeColor','interp','LineWidth',1.4, 'DisplayName','Left outer (mapped)');
        set(ax,'CLim',[1 ncolors]);
        
        % right outer: map to right-inner coords then gradient-plot
        yyaxis(ax,'right'); rightInnerLim = ylim(ax);
        yRouter_mapped = mapLinear(yRouter, rightOuterTrueLim, rightInnerLim);
        h4 = surface(ax, [x;x], [yRouter_mapped;yRouter_mapped], zeros(2,numel(x)), ...
                     repmat(cvals,2,1), ...
                     'FaceColor','none','EdgeColor','interp','LineWidth',1.6, 'DisplayName','Right outer (mapped)');
        set(ax,'CLim',[1 ncolors]);
        
    otherwise
        error('Unknown Mode "%s"', mode);
end

% common axes decorations
xlabel(ax,'$Z_{s^{\prime}}$', 'interpreter', 'latex');
grid(ax,'on'); box(ax,'on');
set(ax,'FontSize', fsize);
view(ax,2); % ensure 2D view for surfaces

% -------------------- overlay axes for outer ticks/labels --------------------
axLout = axes(fig, ...
    'Position', ax.Position, ...
    'Color', 'none', ...
    'XAxisLocation', 'bottom', ...
    'YAxisLocation', 'left', ...
    'Box', 'off', ...
    'XTick', [], ...
    'YColor', [0.6 0 0.6], ... % aesthetic - will be overwritten in perCurve mode below
    'FontSize', fsize, ...
    'HitTest', 'off', ...
    'PickableParts', 'none');
ylabel(axLout, 'Left outer', 'Color', [0.6 0 0.6]);
axLout.TickDir = 'out';
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 30)
set(axLout,'XColor','none');

% Right outer overlay
axRout = axes(fig, ...
    'Position', ax.Position, ...
    'Color', 'none', ...
    'XAxisLocation', 'bottom', ...
    'YAxisLocation', 'right', ...
    'Box', 'off', ...
    'XTick', [], ...
    'YColor', [0 0.6 0.9], ...
    'FontSize', fsize, ...
    'HitTest', 'off', ...
    'PickableParts', 'none');
ylabel(axRout, 'Right outer', 'Color', [0 0.6 0.9]);
axRout.TickDir = 'out';
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 30)
set(axRout,'XColor','none');

% If mode==perCurve, set overlay axis YColor to the assigned outer colours:
if strcmp(mode,'perCurve')
    axLout.YColor = cLouter;
    ylabel(axLout,'$k\frac{d^2\sigma}{dk \, d\Omega_k}$ [barn/sr]', ...
        'Interpreter', 'latex','Color', 'k', 'fontsize', 33);
    axRout.YColor = cRouter;
    ylabel(axRout,'$k\frac{d^2\sigma}{dk \, d\Omega_k}$ [barn/sr]', ...
        'Interpreter', 'latex','Color', 'k', 'fontsize', 33);
else
    % In gradient mode just keep the overlay color neutral; the colormap is visible on lines
end

% position shifts (move tick axes outside)
p = ax.Position;
axLout.Position = [p(1)-leftShift, p(2), p(3), p(4)];
axRout.Position = [p(1)+rightShift, p(2), p(3), p(4)];

uistack(axLout,'top'); uistack(axRout,'top');

% -------------------- update function for zoom/pan/resize --------------------
    function updateOverlays(~,~)
        % keep overlay axes aligned with main axes (handles resize)
        basePos = ax.Position;
        axLout.Position = [basePos(1)-leftShift, basePos(2), basePos(3), basePos(4)];
        axRout.Position = [basePos(1)+rightShift, basePos(2), basePos(3), basePos(4)];
        
        % read current inner Y limits (switch yyaxis reliably)
        yyaxis(ax,'left'); leftInnerLim = ylim(ax);
        nTick = 4;
        yiL = linspace(leftInnerLim(1), leftInnerLim(2), nTick);
        ax.YTick = yiL;
        ax.YTickLabel = arrayfun(@(v)sprintf('%0.f',v), yiL, 'UniformOutput', false);
        yyaxis(ax,'right'); rightInnerLim = ylim(ax);
        yiR = linspace(rightInnerLim(1), rightInnerLim(2), nTick);
        ax.YTick = yiR;
        ax.YTickLabel = arrayfun(@(v)sprintf('%0.f',v), yiR, 'UniformOutput', false);
        
        % remap outer curves into the current inner-limits coordinate system
        switch mode
            case 'perCurve'
                yyaxis(ax,'left');
                yLouter_mapped = mapLinear(yLouter, leftOuterTrueLim, leftInnerLim);
                set(h3, 'YData', yLouter_mapped);
                
                yyaxis(ax,'right');
                yRouter_mapped = mapLinear(yRouter, rightOuterTrueLim, rightInnerLim);
                set(h4, 'YData', yRouter_mapped);
                
            case 'gradient'
                % for surface objects, YData is 2xN, update accordingly
                yyaxis(ax,'left');
                set(h3, 'YData', [ mapLinear(yLouter, leftOuterTrueLim, leftInnerLim); ...
                                   mapLinear(yLouter, leftOuterTrueLim, leftInnerLim) ]);
                % left-inner gradient curve doesn't need remap (native y)
                
                yyaxis(ax,'right');
                set(h4, 'YData', [ mapLinear(yRouter, rightOuterTrueLim, rightInnerLim); ...
                                   mapLinear(yRouter, rightOuterTrueLim, rightInnerLim) ]);
                % right-inner gradient curve doesn't need remap (native y)
        end
        
        % update overlay tick values to reflect TRUE outer scales
        nTick = 5;
        ytL = linspace(leftOuterTrueLim(1), leftOuterTrueLim(2), nTick);
        axLout.YLim = leftOuterTrueLim;
        axLout.YTick = ytL;
        axLout.YTickLabel = arrayfun(@(v)sprintf('%.3g',v), ytL, 'UniformOutput', false);
        
        ytR = linspace(rightOuterTrueLim(1), rightOuterTrueLim(2), nTick);
        axRout.YLim = rightOuterTrueLim;
        axRout.YTick = ytR;
        axRout.YTickLabel = arrayfun(@(v)sprintf('%.3g',v), ytR, 'UniformOutput', false);
    end

% -------------------- attach callbacks so zoom/pan/resize auto-update --------------------
z = zoom(fig); pz = pan(fig);
z.ActionPostCallback = @updateOverlays;
pz.ActionPostCallback = @updateOverlays;
fig.SizeChangedFcn = @updateOverlays;

% also try to add YLim listener (some MATLAB versions allow this)
try
    addlistener(ax, 'YLim', 'PostSet', @updateOverlays);
catch
    % ignore if not supported on older MATLABs
end

% initial sync
updateOverlays();

% -------------------- legend & finishing touches --------------------
% Build a legend using the visible handles (if surfaces used, they are valid handles)
legendHandles = [h1 h2 h3 h4];
Karr=opts.Karr;
legend(legendHandles, {['$k=$', num2str(round(Karr(1)/1e3,3,'significant'))], ...
    ['$k=$', num2str(round(Karr(2)/1e3,3,'significant'))], ...
    ['$k=$', num2str(round(Karr(3)/1e3,3,'significant'))], ...
    ['$k=$', num2str(round(Karr(4)/1e3,3,'significant')), ' [MeV]']}, ...
       'Location','northwest','Orientation','horizontal', 'interpreter', 'latex');
title(ax,'');
set([ax, axLout, axRout], 'FontSize', 25);
grid minor;

end

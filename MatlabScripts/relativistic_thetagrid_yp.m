function [theta,dtheta] = relativistic_thetagrid_yp(Ec,ntheta,display_mode)
%
%	Optimal theta grid for relativistic effects (forward ellipse)
%
%	Input:
%
%		- Ec: kinetic energy of the incoming electron (keV) [1,1] 
%		- ntheta: size of theta grid [1,1] 
%       - display_mode: display grids
%
%	Output: 
%
%		- theta: non-uniform theta grid (rad) [1,ntheta]
%		- dtheta: non-uniform dtheta grid (rad) [1,ntheta]
%
%
%   Dependencies: 
%
%       - pc_dke_yp.m            % physical constants
%       - cumsimps.m             % Cumulative Simpson's method        
%           
% by Y.Peysson CEA-IRFM <yves.peysson@cea.fr>
%
if nargin == 1
    %
    ntheta = 1000;
    display_mode = 0;
    %
elseif nargin == 2
    %
    display_mode = 0;
    %
end
%
[~,~,~,~,~,~,~,mc2,~,~] = pc_dke_yp;%Physics constant
%
gamma = Ec/mc2 + 1;%relativistic gamma
beta = sqrt(1 - 1/gamma^2);%relativistic beta   
%
a = gamma*(1+beta);
b = 1;
%
p = b^2/a;
e = sqrt(a^2 - b^2)/a;
h = b^2/sqrt(a^2 - b^2);
%
hpellipse_ramanujan = pi*(3*(a+b)- sqrt((3*a+b)*(a+3*b)))/2;%Ramanujan formula number 1 for the half ellipse perimeter
%
theta0 = linspace(0,pi,ntheta);
dtheta0 = gradient(theta0);
%
dssdtheta0 = p.*sqrt(1+e^2-2*e*cos(theta0))./(1-e*cos(theta0)).^2;%ellipse elementary arc length
ds0 = dssdtheta0.*dtheta0;
s0 = cumsimps(ds0);
hp_num = s0(end);%numerical half ellipse perimeter
%
if abs(hp_num-hpellipse_ramanujan)/hpellipse_ramanujan > 0.01
    %
    error('Numerical ellipse perimeter is inconsistent with Ramanujan formula number 1. Number of theta points must be increased.')
    %
end
%
s = linspace(0,hp_num,ntheta);
ds = gradient(s);
%
theta = interp1(s0,theta0,s);
dtheta = gradient(theta);
%
theta = theta(:)';
dtheta = dtheta(:)';
%
if display_mode == 1
    %
% Arc length vs theta
figure
plot(theta0/pi, s0, 'rx', ...
     theta/pi, s, 'b*', ...
     theta0/pi, hp_num*ones(size(theta0)), 'k', ...
     'LineWidth', 2, 'MarkerSize', 8);
xlabel('$\theta/\pi$ ', 'Interpreter', 'latex')
ylabel('Arc length [a.u.]', 'Interpreter', 'latex')
legend({'Uniform $\theta$ grid', ...
        'Uniform arc grid', ...
        'Ellipse perimeter (Ramanujan \#1)'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex')
title(['Optimal grid for relativistic effects, $E_c = ', num2str(Ec), '\,\mathrm{keV}$'], ...
      'Interpreter', 'latex')
grid minor
set(gca, 'FontSize', 25)

% dtheta vs theta
figure
plot(theta0/pi, dtheta0, 'k.', 'LineWidth', 0.5, 'MarkerSize', 8, 'displayname', 'Uniform') 
hold on 
plot(theta/pi, dtheta, 'k+','LineWidth', 2, 'MarkerSize', 8, 'displayname', 'Relativistic')
xlabel('$\theta/\pi$ ', 'Interpreter', 'latex')
ylabel('$d\theta$ [rad]', 'Interpreter', 'latex')
title(['Optimal grid for relativistic effects, $E_c = ', num2str(Ec), '\,\mathrm{keV}$'], ...
      'Interpreter', 'latex')
legend('Location', 'NorthWest', 'Interpreter', 'latex')
grid minor
axis([0, 1, 0, 1.1*max(dtheta)])
set(gca, 'YTick', [0.05 0.1], 'YTickLabel', {'0.05','0.1'});
set(gca, 'FontSize', 25)

end

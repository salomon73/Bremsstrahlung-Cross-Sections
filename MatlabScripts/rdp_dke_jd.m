function [seRDP,seRDP_1,seRDP_2] = rdp_dke_jd(ec_in,k_in,t0_in,Z,n2f1,nphi,ntheta,integral_mode,AbsTol,RelTol)
%
%	Doubly differential electron-ion bremstrahlung cross-section (dsigma/dk.domegak)
%	Roche-Ducos-Proriol formula (extended version of Elwer-Haug formula)
%	Ref: E. Haug, Radiation Physics and Chemistry, 207 (2008) 207-214
%
%	Input:
%
%		- ec_in: kinetic energy of the incoming electron (keV) [1,m] 
%		- kc_in: photon energy (keV)  [1,n] 
%		- t0_in: cosine of the angle between the direction of displacement of the incoming  electron and the photon emitted by bremsstrahlung (radian) [1,p] (default = 0)
%		- Z: target ion charge (default = 1) [1,1]
%       - n2f1: number of terms in the series for the calculation o the Gaussian hypergeometric function [1,1]
%         (default : 100)
%       - nphi: number of azimuthal angles for the scattered electron [1,1]
%         (default : 301)
%       - ntheta: number of theta angles for the scattered electron [1,1]
%         (default : 351)
%       - integral_mode: mode of integration over phi and theta of the scatteredd electron ('simps' for Simpson method or 'trapz' or trapeze method)
%         (default : 'simps')
%
%	Output: 
%
%		- seRDP: Full Roche-Ducos-Proriol bremsstrahlung cross-section (m^2) [p,m,n]
%		- seRDP_1: 1st order Roche-Ducos-Proriol bremsstrahlung cross-section, equivalent to Elwert-Haug calculations (m^2) [p,m,n]
%		- seRDP_2: RDP mixed order correction (m^2) [p,m,n]
%		- ec: kinetic energy of the incoming electron (mc2) [p,m,n] 
%		- k: photon energy (mc2)  [p,m,n] 
%		- c: cosine angle between the direction of displacement of the incoming electron and the photon emitted by bremsstrahlung (radian) [p,m,n]
%
%	Warning: Cross-section units : m^2 but energies are in relativistic
%	units.
%
%   Dependencies and toolboxes:
%       
%      Dependencies: 
%
%       - hypergeometric2f1_flint_vec.arch      % MEX compiled hypergeom function (C++)  
%       - cumsimps.m                            % cumulative Simpson's method
%       - pc_dke_yp.m                           % physical constants  
%       - relativistic_thetagrid_yp.m           % relativistic theta grid for polar integration
%       - simps.m                               % Simson's method
%           
%   Toolboxes:
%
%       - MATLAB
%
%
%  by J. Decker EPFL <joan.decker@epfl.ch>, Y.Peysson CEA-DRFC <yves.peysson@cea.fr>
%  and S. Guinchard EPFL <salomon.guinchard@epfl.ch>
%
    if nargin < 4
	    error(2,'Wrong number of input arguments for bhe_dke_yp');
    end
    %
    if nargin == 4
        %
        n2f1 = 100;
        nphi = 101; 
        ntheta = 151; 
        integral_mode = 'simps';
        %
    elseif nargin == 5
        %
        nphi = 101; 
        ntheta = 151; 
        integral_mode = 'simps';
        %
    elseif nargin == 6
        %
        ntheta = 151; 
        integral_mode = 'simps';
        %
    elseif nargin == 7
        %
        integral_mode = 'simps';
        %
    end

    [~,~,~,~,~,~,~,mc2,~,~] = pc_dke_yp; %Physics constant

    if nargin > 8
        if nargin == 9
            RelTol = 1e-6;% default value
        end
        ds3_1 = integral(@(t) ds4(t,nphi,ec_in,k_in,t0_in,Z,n2f1,integral_mode,1),0,pi,'AbsTol',AbsTol,'RelTol',RelTol);%integration over theta
        ds3_2 = integral(@(t) ds4(t,nphi,ec_in,k_in,t0_in,Z,n2f1,integral_mode,2),0,pi,'AbsTol',AbsTol,'RelTol',RelTol);%integration over theta
    else

        phi = linspace(0,2*pi,nphi);
        theta = relativistic_thetagrid_yp(max(ec_in),ntheta);

        if exist('integral_mode', 'var') && strcmp(integral_mode,'simps') && exist('simps')
            %
            ds3_1 = simps(theta,ds4(theta,phi,ec_in,k_in,t0_in,Z,n2f1,integral_mode,1),2);%integration over theta
            ds3_2 = simps(theta,ds4(theta,phi,ec_in,k_in,t0_in,Z,n2f1,integral_mode,2),2);%integration over theta
            %
        else
            %
            ds3_1 = trapz(theta,ds4(theta,phi,ec_in,k_in,t0_in,Z,n2f1,integral_mode,1),2);%integration over theta
            ds3_2 = trapz(theta,ds4(theta,phi,ec_in,k_in,t0_in,Z,n2f1,integral_mode,2),2);%integration over theta
            %
        end

    end
    %
    mask = ec_in > k_in; % Photon are only emitted by electron of higher energies
    %
    seRDP_1 = ds3_1.*mask;
    seRDP_2 = ds3_2.*mask;
    %
    test = abs(squeeze(imag(seRDP_1)).*k_in/Z^2/mc2/1e-31);% millibarns (10^-31 m2)
    %
    if find(test > eps) 
        %
        warning('Some values of the first order RDP cross-section have an imaginary part larger than numerical accuracy.')
        %
    end
    %
    seRDP = real(seRDP_1) + seRDP_2;
    %
    seRDP(isnan(seRDP)) = 0;
    %
end
    
function ds4 = ds4(theta_in,np,ec_in,k_in,c0_in,Z,n2f1,integral_mode,order)
    %
    if isnan(np)
        ds4 = integral(@(phi) ds5(theta_in,phi,ec_in,k_in,c0_in,Z,n2f1,integral_mode,order),0,2*pi,'RelTol',1e-7,'ArrayValued',true);%integration over phi
    else
        phi_in = linspace(0,2*pi,np);
        %
        if exist('integral_mode', 'var') && strcmp(integral_mode,'simps') && exist('simps')
            %
            ds4 = simps(phi_in,ds5(theta_in,phi_in,ec_in,k_in,c0_in,Z,n2f1,integral_mode,order),1);%integration over phi
            %
        else
            %
            ds4 = trapz(phi_in,ds5(theta_in,phi_in,ec_in,k_in,c0_in,Z,n2f1,integral_mode,order),1);%integration over phi
            %
        end
    end
end

function ds5 = ds5(theta_in,phi_in,ec_in,k_in,c0_in,Z,n2f1,integral_mode,order)    
    %
    [~,~,~,~,~,~,re,mc2,~,alpha] = pc_dke_yp; %Physics constant
    %
    nt = length(theta_in);
    np = length(phi_in);
    %
    pones = ones(np,1);
    tones = ones(1,nt);
    ptones = ones(np,nt);
    %
    phi = phi_in(:)*tones;
    %
    ec = ec_in/mc2*ptones; % Relativistic units
    k = k_in/mc2*ptones;   % Relativistic units
    c0 = c0_in*ptones;
    c = pones*cos(theta_in); 
    %
    s0 = sqrt(1 - c0.^2);
    s  = sqrt(1 - c.^2);
    %
    e0 = ec + 1;          % Initial total electron energy
    p0 = sqrt(e0.^2 - 1); % Initial electron momentum
    ep = e0 - k;          % Final electron energy
    p  = sqrt(ep.^2 - 1); % Final electron momentum
    %
    q2 = p.^2 + p0.^2 + k.^2 - 2*p0.*k.*c0 + 2*p.*k.*c - 2*p0.*p.*(c.*c0 + s.*s0.*cos(phi));
    %
    a0 = alpha*Z*e0./p0;
    a  = alpha*Z*ep./p;
    %
    N = 4*pi^2*a0.*a./(exp(2*pi*a0)-1)./(1-exp(-2*pi*a));
    %
    D0 = 2*(e0.*k-p0.*c0.*k);
    D  = 2*(ep.*k-p.*c.*k);
    %
    mhu = 2*(e0.*ep+p0.*p-1);
    %
    x = 1 - D0.*D./mhu./q2;
    %
    PHI = a0.*log(q2./D) - a.*log(mhu./D);

    % matlab built in hypergeom calls
    if ~isnan(n2f1)
        V = gaussianhypergeometric(-1i*a0,1i*a,1,x,n2f1);
        W = gaussianhypergeometric(1-1i*a0,1+1i*a,2,x,n2f1);
    else
        V = reshape(hypergeometric2f1_flint_vec(-1i*a0, 1i*a, 1*ones(size(x)), x),size(x));
        W = reshape(hypergeometric2f1_flint_vec(1-1i*a0, 1+1i*a, 2*ones(size(x)), x),size(x));
    end
    %
    VWsq2  = (V+1i*a.*x.*W)./q2;
    xWsDD0 = (1-x).*W./(D0.*D);
    mD0p1  = mhu./D0+1;
    mDm1   = mhu./D-1;
    %
    J1 = 2*(VWsq2.*(ep./D0 - e0./D) + 1i.*xWsDD0.*(e0.*a.*mDm1 - ep.*a0.*mD0p1));
    J1c = conj(J1);
    %
    Xk  = {zeros(size(k)),zeros(size(k)),k};
    Xkn = {zeros(size(k)),zeros(size(k)),ones(size(k))};
    %
    Xp0 = {p0.*s0,zeros(size(p0)),p0.*c0};
    %
    Xp = {p.*s.*cos(phi),p.*s.*sin(phi),p.*c};
    %
    Xq = {Xp0{1} - Xp{1} - Xk{1},Xp0{2} - Xp{2} - Xk{2},Xp0{3} - Xp{3} - Xk{3}};
    %
    XP = {p0.*Xp{1} + p.*Xp0{1},p0.*Xp{2} + p.*Xp0{2},p0.*Xp{3} + p.*Xp0{3}};
    %
    XJ2 = {(VWsq2./D).*Xq{1} - 1i*a.*xWsDD0.*(mDm1.*Xq{1} - XP{1}./p0),...
           (VWsq2./D).*Xq{2} - 1i*a.*xWsDD0.*(mDm1.*Xq{2} - XP{2}./p0),...
           (VWsq2./D).*Xq{3} - 1i*a.*xWsDD0.*(mDm1.*Xq{3} - XP{3}./p0)};
    %
    XJ2c = {conj(XJ2{1}), conj(XJ2{2}), conj(XJ2{3})};

    % man16
    % XJ3 = {(VWsq2./D).*Xq{1} - 1i*a0.*xWsDD0.*(mD0p1.*Xq{1} - XP{1}./p),...
    %        (VWsq2./D).*Xq{2} - 1i*a0.*xWsDD0.*(mD0p1.*Xq{2} - XP{2}./p),...
    %        (VWsq2./D).*Xq{3} - 1i*a0.*xWsDD0.*(mD0p1.*Xq{3} - XP{3}./p)};

    % man17 hau08
    XJ3 = {(VWsq2./D0).*Xq{1} - 1i*a0.*xWsDD0.*(mD0p1.*Xq{1} - XP{1}./p),...
           (VWsq2./D0).*Xq{2} - 1i*a0.*xWsDD0.*(mD0p1.*Xq{2} - XP{2}./p),...
           (VWsq2./D0).*Xq{3} - 1i*a0.*xWsDD0.*(mD0p1.*Xq{3} - XP{3}./p)};
    %
    XJ3c = {conj(XJ3{1}), conj(XJ3{2}), conj(XJ3{3})};

    fopt=1; % (1) hau08 (2) man17 man16 (3) elw69
    
    if order == 1
        %
        A1 = (e0.*ep - 1 - sdot(Xkn,Xp0).*sdot(Xkn,Xp)).*J1.*J1c;
        A2 = (e0.*ep + 1 + sdot(Xkn,Xp0).*sdot(Xkn,Xp)).*(sdot(XJ2,XJ2c) + sdot(XJ3,XJ3c));
        if fopt == 1 || fopt == 3
            A3 =  2*real((sdot(XJ3c,Xp0) - sdot(XJ2c,Xp0)).*sdot(XJ2,Xkn).*sdot(Xp,Xkn) - (sdot(XJ3c,Xp) - sdot(XJ2c,Xp)).*sdot(XJ3,Xkn).*sdot(Xp0,Xkn));
        elseif fopt == 2
            A3 =  2*real((sdot(XJ3,Xp0) - sdot(XJ2,Xp0)).*sdot(XJ2,Xkn).*sdot(Xp,Xkn) - (sdot(XJ3,Xp) - sdot(XJ2,Xp)).*sdot(XJ3,Xkn).*sdot(Xp0,Xkn));
        end
        if fopt == 1 || fopt == 2
            A4 = -2*real((e0.*ep + 1 + sdot(Xp0,Xp)).*sdot(XJ2,Xkn).*sdot(XJ3c,Xkn));
            A5 =  2*real(sdot(XJ2,Xp0).*sdot(XJ3c,Xp) - sdot(XJ2,Xp).*sdot(XJ3c,Xp0));
        elseif fopt == 3
            A4 = -2*real((e0.*ep + 1 + sdot(Xp0,Xp)).*sdot(XJ2c,Xkn).*sdot(XJ3,Xkn)); 
            A5 =  2*real(sdot(XJ2c,Xp0).*sdot(XJ3,Xp) - sdot(XJ2,Xp).*sdot(XJ3c,Xp0));
        end   
        A6 =  2*real(e0.*J1c.*(sdot(XJ3,Xp) - sdot(XJ2,Xkn).*sdot(Xp,Xkn)));
        A7 =  2*real(ep.*J1c.*(sdot(XJ2,Xp0) - sdot(XJ3,Xkn).*sdot(Xp0,Xkn)));
        if fopt == 1
            A8 =  2*real(sdot(XJ2,XJ3c).*(sdot(Xp0,Xp) - sdot(Xkn,Xp0).*sdot(Xkn,Xp)));
        elseif fopt == 2
            A8 =  2*real(sdot(XJ2,XJ3).*(sdot(Xp0,Xp) - sdot(Xkn,Xp0).*sdot(Xkn,Xp)));
        elseif fopt == 3
            A8 =  2*real(sdot(XJ2c,XJ3).*(sdot(Xp0,Xp) - sdot(Xkn,Xp0).*sdot(Xkn,Xp)));
        end
        %
        d5 = alpha*Z^2*re^2*k.*p.*N.*(A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8)./p0/pi^2;
    elseif order == 2
        %
        B11 = (e0.*ep - 1 - sdot(Xkn,Xp0).*sdot(Xkn,Xp)).*J1;
        B12 =  e0.*(sdot(XJ3,Xp) - sdot(XJ2,Xkn).*sdot(Xp,Xkn));
        B13 =  ep.*(sdot(XJ2,Xp0) - sdot(XJ3,Xkn).*sdot(Xp0,Xkn));
        % 
        B1 = 2*real((B11 + B12 + B13).*(cos(PHI) + 1i*sin(PHI)))/pi;
        B2 = alpha*Z*sdot(Xq,Xk).*(e0.*ep - 1 - sdot(Xkn,Xp0).*sdot(Xkn,Xp))./sqrt(q2)./D0./D;
        %
        d5 = alpha^2*Z^3*re^2*k.*p.*N.*sdot(Xq,Xk).*(B1 + B2)./p0./sqrt(q2)./D0./D;
    end
    %
    ds5 = d5.*s; % sin theta
end
%
function [y] = gaussianhypergeometric(a,b,c,z,n2f1)
    %
    % Resampling for vectorial calculations or by mex function (specific for the bremsstrahlung cross-section calculation)
    %
    y = hypergeometric2f1(a(:),b(:),c*ones(length(z(:)),1),z(:),n2f1);
    %
    y = reshape(y,size(a));
    %
end
%
function [s] = sdot(X,Y)
    %
    % generalized scalar product for vectorial calculations
    %
    s = X{1}.*Y{1} + X{2}.*Y{2} + X{3}.*Y{3};
    %
end


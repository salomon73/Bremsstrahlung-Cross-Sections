function [seBHE,seBH,Elwert,ec,k,c] = bhe_dke_yp(ec_in,k_in,t0_in,Z,params,smodel)
%
%	This script calculates the unscreened DDCS in photon energy and solid
%	angle from the Bethe-Heitler TDCS 
%   (1BS, 2BN formulas in H.W.Koch and J.W.Motz, Rev. Mod. Phys. 31, 4 (1959) 920) 
%   including Elwert factor (Coulomb), and screening factor F(q).  
%
%	Input:
%
%		- ec_in: kinetic energy of the incoming electron (keV) [1,m] 
%		- kc_in: photon energy (keV)  [1,n] 
%		- t0_in: cosine of the angle between the direction of displacement of the incoming 
%		  electron and the photon emitted by bremsstrahlung (radian) [1,p]
%		  (default = 0)
%		- Z: target ion charge (default = 1) [1,1]
%       - params = calculation parameters structure (for 1BS -> 2BN integration) 
%       - smodel = screening model structure
%           (default = [])
%
%	Output: 
%
%		- seBHE: Bethe-Heitler + Elwert factor bremstrahlung cross-section (m^2) [p,m,n]
%		- seBH: Bethe-Heitler bremstrahlung cross-section (m^2) [p,m,n]
%		- Elwert: Elwert correction factor [p,m,n]
%		- ec: kinetic energy of the incoming electron (mc2) [p,m,n] 
%		- k: photon energy (mc2)  [p,m,n] 
%		- c: cosine angle between the direction of displacement of the incoming electron and the photon emitted by bremsstrahlung (radian) [p,m,n]
%
%	Warning: Cross-section units : m^2 but energies are in relativistic units. To get 
%			 cross-sections in standard m^2/keV units, seBH or seBHe must be divided 
%			 by 511 keV.
%
%   Dependencies: 
%
%       - cumsimps.m                      % Cumulative Simpson      
%       - fionradius_yp.m                 % Effective ion radius         
%       - load_yukawa.m                   % Multi-Yukawa model parameters      
%       - pc_dke_yp.m                     % Physical constants      
%       - relativistic_thetagrid_yp.m     % Relat' grid for theta                      
%       - simps.m                         % Simpson's method 
%       - dft_integrals_yp.m              % Density Functional Theory              
%       - tf_potential_yp.m               % Thomas-Fermi atomic model           
%       - trapz_dke_yp.m                  % Trapeze method   
%       - MYi_Zj.m                        % DFT multi-Yukawa coefficients       
%
% by Y.Peysson CEA-DRFC <yves.peysson@cea.fr>, J. Decker MIT-RLE <jodecker@mit.edu>
% and S. Guinchard EPFL <salomon.guinchard@epfl.ch>
%
if nargin < 4
	error(2,'Wrong number of input arguments for bhe_dke_yp');
end
%
if nargin < 5
    params = struct;
	smodel =  struct; %no screening
end
%
if nargin < 6
	smodel =  struct; %no screening
end
%
if isempty(smodel) && isempty(params)
    params = struct;
	smodel =  struct; %no screening
end
%
[~,~,~,~,~,~,re,mc2,~,alpha] = pc_dke_yp; %Physics constant
%
if isempty(fields(smodel)) || isempty(fields(params))
    %
    % Calculation of 2BN formula in H.W.Koch and J.W.Motz, Rev. Mod. Phys. 31, 4 (1959) 920
    % Bethe-Heitler DDCS (analytical)
    %
    ec = repmat(ones(length(t0_in),1)*ec_in(:)'/mc2,[1,1,length(k_in)]);            %Relativistic units
    k = shiftdim(repmat(ones(length(ec_in),1)*k_in(:)'/mc2,[1,1,length(t0_in)]),2); %Relativistic units
    if (isscalar(ec_in)) && (isscalar(t0_in)) % scalar case
        k = reshape(k,1,1,length(k_in));
    end
    c = repmat(t0_in(:)*ones(1,length(ec_in)),[1,1,length(k_in)]);
    %
    s = sqrt(1-c.^2);
    %
    mask = ec > k;  %Photon are only emitted by electron of higher energies
    %
    e0 = ec+1;
    p0 = sqrt(e0.^2-1);
    p = sqrt((e0-k).^2-1);
    ep = e0-k;d0 = e0-c.*p0;
    %
    ksi0 = alpha*Z*e0./p0;
    ksi1 = alpha*Z*ep./p;
    %
    Elwert = (ksi1./ksi0).*((1-exp(-2*pi*ksi0))./(1-exp(-2*pi*ksi1)));
    %
    q   = sqrt(p0.^2+k.^2-2*p0.*k.*c);
    w1  = log((ep.*e0-1+p.*p0)./(ep.*e0-1-p.*p0));
    w2  = log((ep+p)./(ep-p));
    w3  = log((q+p)./(q-p));
    s1  = 8*(s.^2).*(2*e0.^2+1).*(p0.^(-2)).*(d0.^(-4));
    s2  = 2*(5*e0.^2+2*e0.*ep+3).*(p0.^(-2)).*(d0.^(-2));
    s3  = 2*(p0.^2-k.^2).*(q.^(-2)).*(d0.^(-2));
    s4  = 4*ep.*(p0.^(-2)).*(d0.^(-1));
    s5  = 4*e0.*(s.^2).*(3*k-(p0.^2).*ep).*(p0.^(-2)).*(d0.^(-4));
    s6  = 4*(e0.^2).*(e0.^2+ep.^2).*(p0.^(-2)).*(d0.^(-2));
    s7  = (2-2*(7*e0.^2-3*e0.*ep+ep.^2)).*(p0.^(-2)).*(d0.^(-2));
    s8  = 2*k.*(e0.^2+ep.*e0-1).*(p0.^(-2)).*(d0.^(-1));
    s9  = 4*w2.*(p.^(-1)).*(d0.^(-1));
    s10 = w3.*(p.^(-1)).*(q.^(-1));
    s11 = 4*d0.^(-2);
    s12 = 6*k.*d0.^(-1);
    s13 = 2*k.*(p0.^2-k.^2).*(q.^(-2)).*(d0.^(-1));
    s14 = w1.*(p.^(-1)).*(p0.^(-1));
    %
    s0 = s1-s2-s3+s4+s14.*(s5+s6+s7+s8)-s9+s10.*(s11-s12-s13);
    a = alpha*Z^2*re^2*p./(8*pi*k.*p0); %Cross-section in m^2 unit
    %
    seBH  = s0.*a.*mask;
    seBHE = seBH.*Elwert;
    %
    seBHE(isnan(seBHE))   = 0;
    seBHE(~isreal(seBHE)) = 0;
else
    %
    % With Atomic Screening specified by smodel.screening_model
    %
    if strcmp(smodel.screening_model,'moliere')
        % Moliere's screening model (fully analytical)
        % based on the Thomas-Fermi atomic model
        % TODO
        [seBHE,seBH,Elwert,ec,k,c] = bhe_moliere(ec_in,k_in,t0_in,Z); 
        %
    else
        % All other screening models (~Molière)
        [seBHE,seBH,Elwert,ec,k,c] = bhe_dke_kk(ec_in,k_in,t0_in,Z,params,smodel); 
    end
end

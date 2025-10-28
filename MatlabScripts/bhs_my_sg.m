function bhs_my_sg = bhs_my_sg(E0_in, k_in, theta0, Z0, YKorder, nion, A, b_bar, c_bar2)
%
% This function computes the screened Bethe-Heitler bremsstrahlung cross section for 
% electrons hitting a target of atomic number Z0, for arbitrary ionisation states 
% The screening is taken into account with a Multi-Yukawa atomic
% model whose order is specified in input. 
% 
%
%   Input: 
%
%       - E0_in   : Energy of incoming electron (keV) 
%       - k_in    : array of energies of emitted photons (keV)
%       - theta0  : angle between incoming e and photon (rad)
%       - Z0      : atomic number of the target element
%       - params  : struct with a field nion for a slice of ionization
%                    states to compute.
%       - YKorder : order of the multi-Yukawa model.
%       - nion    : slice of ionized states
%       - A       : array of Yukawa weights
%       - b_bar   : array of inverse screening lengths
%       - c_bar2  : array of inv screening lengths * ionization degree qZ
%
%   Output:
%   
%       - bhs_my_sigma : doubly differential cross section integrated 
%                        over electron degrees of freedom (theta, phi)
%
%   Dependencies and toolboxes:
%       
%      Dependencies: 
%
%       - pc_dke_yp.m            % physical constants
%       - MYi_Zj.m               % DFT multi-Yukawa coefficients   
%       - bhe_dke_yp.m           % For Bethe--Heitler--Sauter      
%           
%   Toolboxes:
%
%       - MATLAB
%
%   By S. Guinchard <salomon.guinchard@epfl.ch> 
%   Last update (28/10/2025)



% Define necessary quantities
%_____________________________

[~,~,~,~,~,~,re,mc2,~,alpha] = pc_dke_yp; % Physics constant

% Normalize shapes: ensure k_in is row vector
k_in = k_in(:)';    % 1 x nint

% relativistic units (scalar E0, vector k)
E0 = E0_in./mc2 + 1;        % scalar (total energy / mc^2)
p0 = sqrt(E0.^2 - 1);       % scalar
k_in = k_in./mc2;           % 1 x nint
E = E0 - k_in;              % 1 x nint
p = sqrt(E.^2 - 1);         % 1 x nint


% geometry (scalar theta0)
c0 = cos(theta0);
s0 = sin(theta0);

% p0·k = p0 * k * cos(theta0) (1 x nint)
dotkp0 = p0 .* k_in .* c0;                    % 1 x nint

% |p0 x k|^2 = (p0*k*sin theta0)^2
modp0crossk2 = (p0 .* k_in .* s0).^2;         % 1 x nint


% p0 vector = magnitude p0 along direction (cos theta0, 0, sin theta0)
% |p0_vec - k_vec| assuming k along z axis 
modp0mk = sqrt( (p0.*c0 - k_in).^2 + (p0.*s0).^2 ); % 1 x nint

D1 = 2*(E0.*k_in - dotkp0);

% Compute BH DDCS 
%________________

fact = alpha*Z0^2*re^2/2/pi./k_in./p0;
[~,seBH(1,1,:),~,~,~,~] = bhe_dke_yp(E0_in,mc2.*k_in,cos(theta0), Z0);   
seBH = squeeze(seBH)';

% Compute DDCS 
%_____________
nint = numel(k_in);
sigma = zeros(numel(nion), nint);

for kk = 1:numel(nion)
    for i = 1:YKorder

        % Diagonal term
        b = b_bar(kk,i);
        diag_term = A(kk,i)^2 * ( ...
                c_bar2(kk,i)^2 * (I20(E0,E,k_in,p,D1,modp0crossk2,modp0mk,b,seBH, fact)) + ...
            2 * c_bar2(kk,i)   * (I21(E0,E,k_in,p,D1,modp0crossk2,modp0mk,b, fact)) + ...
                                 I1(E0,E,k_in,p,D1,modp0crossk2,modp0mk,b, fact));
        sigma(kk,:) = sigma(kk,:) + diag_term;

        % Non-diagonal terms
            for j = i+1:YKorder
                
                bj = b_bar(kk,j); bi = b_bar(kk,i);
                denom = bj^2 - bi^2;

                factor = A(kk,i)*A(kk,j) / denom;

                cross_sum = c_bar2(kk,i)*c_bar2(kk,j)*( ...
                    I10(E0,E,k_in,p,D1,modp0crossk2,modp0mk,bi,seBH, fact) - I10(E0,E,k_in,p,D1,modp0crossk2,modp0mk,bj,seBH, fact)) + ...
                    ...
                    (c_bar2(kk,i)+c_bar2(kk,j))*( ...
                    I11(E0,E,k_in,p,D1,modp0crossk2,modp0mk,bi, fact) - I11(E0,E,k_in,p,D1,modp0crossk2,modp0mk,bj, fact)) + ...
                    ...
                    (I2(E0,E,k_in,p,D1,modp0crossk2,modp0mk,bi, fact) - I2(E0,E,k_in,p,D1,modp0crossk2,modp0mk,bj,fact));

                non_diag_term = 2 * factor * cross_sum;
                sigma(kk,:) = sigma(kk,:) + non_diag_term;
            end
    end
end
% final DDCS
bhs_my_sg = sigma;
end

% Fronsdal-Uberall-Haug 1st integral
%___________________________________

function I1 = I1(E0,E,k,p,D1,modp0crossk2,modp0mk,b, fact)

    W  = sqrt((E0.*D1./k - 2 + b^2).^2 + 4.*modp0crossk2./k.^2);
    L1 = log(((E0.*E - 1).*D1./k + E*b^2 + p.*W)./((E0.*E - 1).*D1./k + E*b^2 - p.*W));
    L2 = log(((modp0mk + p).^2 + b^2)./((modp0mk - p).^2 + b^2));
    %
    I1_1 = 16.*p.*(4.*E0.^2 + b^2).*modp0crossk2./k.^2./W.^4;
    %
    I1_2 = -(2*p.*(4*E0.^2 + 2*E0.*E - E.*D1./k + b^2*(1 - 2*E.*k./D1))./W.^2);
    %
    I1_3 = 2*p.*( ...
             ( D1.^2 + 2*(E0.*E -1).*D1 + b^2*(D1 - 2*E.*k))./((D1 + b^2).^2 + 4*p.^2.*b^2) ...
          .* ( (16*E0.*E - 4*E0.^2*b^2 - b^4)./D1 - (4*E0.^2 + b^2).*(D1.^2 + 2*(E0.*E - 1).*D1 + b^2*(D1 - 2*E.*k))./k.^2./W.^2) ...
        )./W.^2;
    %
    I1_4 = -(2*k.^2.*p).*(4*(4*E.^2+(1-D1).*b.^2)./D1.^2 + (D1.^2 + 2*(E0.*E - 1).*D1 + b^2*(D1 - 2*E.*k))./D1./modp0mk.^2)./((D1 + b^2).^2 + 4*p.^2.*b^2);
    %
    I1_5 = -4*k.*log(E+p)./D1;
    %
    I1_6 = L1.*(...
           2*k + 4*k.*(E0.^2 + p.^2 + b^2)./D1 ...
        +  2*(E0.*D1 - 2*k + b^2*k).*(8*E0.*E - D1.^2/2 - b^2*(2*E0.^2+ 2*p.^2 + D1) - b^4)./D1./W.^2 ...
        + (2*(2*E0.^2 + b^2).*(D1 - 2*E.*k) + D1.^2 + 2*(E0.*E - 1).*D1)./k./W.^2 ...
        -  3*(4*E0.^2 + b^2).*(E0.*D1./k - 2 + b^2).*(D1.^2 + 2*(E0.*E - 1).*D1 + b^2.*(D1 - 2*E.*k))./k./W.^4 ...
        )./W;

    I1_7 = k.^2.*L2.*(2./D1 - 2 + (D1 -2*E.*k)./2./modp0mk.^2)./D1./modp0mk;
    %
    I1 = fact .* (I1_1 + I1_2 + I1_3 + I1_4 + I1_5 + I1_6 + I1_7);

end

% Fronsdal-Uberall-Haug 2nd integral
%___________________________________

function I2 = I2(E0,E,k,p,D1,modp0crossk2,modp0mk,b, fact)
    %
    W  = sqrt((E0.*D1./k - 2 + b^2).^2 + 4.*modp0crossk2./k.^2);
    L1 = log(((E0.*E - 1).*D1./k + E*b^2 + p.*W)./((E0.*E - 1).*D1./k + E*b^2 - p.*W));
    L2 = log(((modp0mk + p).^2 + b^2)./((modp0mk - p).^2 + b^2));
    %
    I2_1 = 2*p.*(4*E0.^2 + b^2).*(E0.*D1./k - 2 + b^2)./W.^2;
    I2_2 = L1.*(k.*(D1+2*b^2) + 2*k.*(b^4 + 2*b^2*(E0.^2+p.^2) - 8*E0.*E)./D1 + (4*E0.^2 + b^2).*((D1 + 2*E0.*E -2).*D1 + b^2*(D1 - 2*E.*k))./k./W.^2)./W;
    I2_3 = k.^2.*L2.*(2*(4*E.^2 + b^2*(1-D1))./D1 + (D1.^2 + 2*D1.*(E0.*E - 1) +b^2*(D1 - 2* E.*k))./modp0mk.^2/2)./D1./modp0mk;
    I2_4 = 4*k*b^2.*log(E+p)./D1;
    %
    I2 = - fact.*( I2_1 + I2_2 + I2_3 - I2_4);
end

% Ionization terms
%_________________

function I10 = I10(E0,E,k,p,D1,modp0crossk2,modp0mk,b, seBH, fact)
    I10 = (I2(E0,E,k,p,D1,modp0crossk2,modp0mk,b, fact) - I2(E0,E,k,p,D1,modp0crossk2,modp0mk,0, fact) ...
           + b^2.*seBH)./b.^4;
end

function I20 = I20(E0,E,k,p,D1,modp0crossk2,modp0mk,b, seBH, fact)
    I20 = 2*(I2(E0,E,k,p,D1,modp0crossk2,modp0mk,b, fact) - I2(E0,E,k,p,D1,modp0crossk2,modp0mk,0, fact)+b^2.*seBH)./b.^6 ...
        + (I1(E0,E,k,p,D1,modp0crossk2,modp0mk,b, fact)-seBH)./b.^4;
end

function I11 = I11(E0,E,k,p,D1,modp0crossk2,modp0mk,b, fact)
    I11 = (I2(E0,E,k,p,D1,modp0crossk2,modp0mk,0, fact) - I2(E0,E,k,p,D1,modp0crossk2,modp0mk,b, fact))./b.^2;
end

function I21 = I21(E0,E,k,p,D1,modp0crossk2,modp0mk,b, fact)
    I21 = - I1(E0,E,k,p,D1,modp0crossk2,modp0mk,b, fact)./b.^2 ...
          + (I2(E0,E,k,p,D1,modp0crossk2,modp0mk,0, fact) - I2(E0,E,k,p,D1,modp0crossk2,modp0mk,b, fact))./b.^4;
end

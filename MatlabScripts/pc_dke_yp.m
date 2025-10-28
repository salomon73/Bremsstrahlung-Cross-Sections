function [varargout] = pc_dke_yp(~)
%
% Calculates the universal physics constants
% (see www.tcaep.co.uk/science/constant)
%
%   INPUTS: None
%
%
%   OUTPUTS: [qe,me,mp,mn,e0,mu0,re,mc2,clum,alpha,kB]
%   
%   - qe : Absolute value of the electron charge (C)
%   - me : Electron rest mass (Kg)
%   - mp : Proton rest mass (Kg)
%   - mn : Neutron rest mass (Kg)
%   - e0 : Free space permeability (F/m)
%   - mu0 : Free space magnetic permeability (H/m)
%   - re : Classical electron radius (m)
%   - mc2 : Electron energy rest mass (keV)
%   - clum : Speed of light (m/s)
%   - alpha : Fine structure constant
%   - kB : Boltzmann's constant (J/K-1)
%   - amu : Atomic mass unit (Kg)
%   - h : Planck constant
%   - hbar : Reduced Planck constant
%   - a0 : Bohr radius (m)
%   - lambdaC : Compton length (m)
%
% by Joan Decker <jodecker@mit.edu> (MIT/RLE) and Yves Peysson <yves.peysson@cea.fr> (CEA/DRFC)
%
qe = 1.602176462e-19;%Absolute value of the electron charge (C)
me = 9.10938188e-31;%Electron rest mass (Kg)
mp = 1.67262158e-27;%Proton rest mass (Kg)
mn = 1.67492716e-27;%Neutron rest mass (Kg)
e0 = 8.854187818e-12;%Free space permittivity (F/m)
mu0 = 12.566370614e-7;%Free space permeability (F/m)
re = 2.817940285e-15;%Classical electron radius (m)
mc2 = 510.998902;%Electron energy rest mass (keV)
clum = 299792458;%Speed of light (m/s)
alpha = 7.297352533e-3;%Fine structure constant (~1/137)
kB = 1.3806504e-23;%Boltzmann's constant (J/K-1)
amu = 1.66053906660e-27;%Atomic mass unit (Kg)
h =  6.62607015e-34;% J.s
%
if nargin > 0,%CGS units
    qe = qe*clum*1e1;% statcoul
    me = me*1e3;% g
    mp = mp*1e3;% g
    mn = mn*1e3;% g
    e0 = NaN;% does not exist in CGS units
    mu0 = NaN;% does not exist in CGS units
    re = re*1e2;% cm
    clum = clum*1e2;% cm/s
    mc2 = mc2*qe*1e11/clum;% erg
    alpha = alpha;
    kB = kB*1e7; % erg/K
    amu = amu*1e3;% g
    h =  6.62607015e-27;% erg.s
end
%
hbar = h/2/pi;%(J.s or Erg.s)
lambdaC = h/me/clum;%Compton wavelength (m or cm)
a0 = lambdaC/2/pi/alpha;%Bohr radius (m or cm)
%
if nargout >= 1, varargout{1} = qe;end
if nargout >= 2, varargout{2} = me;end
if nargout >= 3, varargout{3} = mp;end
if nargout >= 4, varargout{4} = mn;end
if nargout >= 5, varargout{5} = e0;end
if nargout >= 6, varargout{6} = mu0;end
if nargout >= 7, varargout{7} = re;end
if nargout >= 8, varargout{8} = mc2;end
if nargout >= 9, varargout{9} = clum;end
if nargout >= 10, varargout{10} = alpha;end
if nargout >= 11, varargout{11} = kB;end
if nargout >= 12, varargout{12} = amu;end
if nargout >= 13, varargout{13} = h;end
if nargout >= 14, varargout{14} = hbar;end
if nargout >= 15, varargout{15} = lambdaC;end
if nargout >= 16, varargout{16} = a0;end


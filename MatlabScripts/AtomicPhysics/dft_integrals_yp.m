function [a_DFT,out] = dft_integrals_yp(rn,dens,Z,N,Z0,model)
%
% Calculation of the integrals I1 and I2 for the Fokker-Planck screening function g, and the determination of the effective length (asymptotic high momentum
% limit) from Thomas-Fermi or Hartree-Fock-Slater models
%
%   INPUTS:
%
%       - rn : normamized atom radius (Bohr or Thomas-Fermi radius)
%       - dens : normalized electron density units
%       - Z : Atomic number of the considered element 
%       - N : Number of bounded electrons 
%       - Z0 : charge of the ionized atom Z0 = Z - N
%       - model : DFT model ('Thomas-Fermi','Tseng-Pratt')
%
%   OUTPUTS:
%
%       - a_DFT : output structure containing the ion radius (bohr radius units), model, and electron distribution
%       - out : output structure of the numerical solution (out.I1 and out.I2 Hesslow's integrals and all ion radii, and Thomas-Fermi and Tseng-Pratt radii)
%
%   By Y. Peysson CEA-IRFM <yves.peysson@cea.fr>, A. Jardin IJF/PAN <axel.jardin@ifj.edu.pl>, Krzysztof Krol IFJ/PAN <krzysztof.krol@ifj.edu.pl>
%
if nargin < 5
    %
    error('Not enough input arguments in dft_integrals_yp');
    %
end
%
if nargin < 6
    %
    model = 'Thomas-Fermi';
    %
end
%
[~,~,~,~,~,~,~,~,~,alpha,~,~] = pc_dke_yp;%Universal physics constants
%
gammaE = 0.5772156649;%Euler-Mascheroni constant
%
rn = rn(:);
dens = dens(:);
dens(1)=dens(2);%avoid singularity in DFT calculation for the first point
%
out.N = 4*pi*trapz_dke_yp([rn,dens.*rn.^2]);%Should be equal to N
%
err_N = abs(out.N - N)/N;
%
if err_N > 0.05
    %
    disp(['WARNING: Numerical accuracy becomes poor for the number of bounded N electrons: ',num2str(out.N),' for ',num2str(N),'.']);
    %
end
%
if rn(1) == 0, rn(1) = eps;end
%
out.I1 = 4*pi*simps(rn,dens.*rn.^2.*log(rn))/out.N;%formula B5 (Hesslow's phd thesis) -> formula B4
%
for is = 1:length(rn)
    %
    s = rn(is);
    %t = linspace(eps,s,1000);
    t = rn(1:is);
    %
    dens_p = interp1(rn,dens,(s+t)/2);
    %
    dens_m = interp1(rn,dens,(s-t)/2);
    dens_m(isnan(dens_m)) = dens(1);
    kernel_t = (s.^2 - t.^2).*(s.^2.*log(s) - t.^2.*log(t)).*dens_p.*dens_m;
    %
    if length(kernel_t) ==  1
        %
        kernel_s(is) = 0;
        %
    else
        %
        kernel_s(is) = simps(t(:),kernel_t(:));
        %
    end
    %
end
%
kernel_s(isnan(kernel_s)) = 0;
%
out.I2 = (pi/out.N)^2*simps(rn(:),kernel_s(:));%formula B8 (Hesslow's phd thesis) -> formula B7
%
out.a_TF_DFT = 2*exp(gammaE - 1+(2*Z*out.I1 + out.N*(7/6 - out.I2))/(Z+Z0))/alpha;%Thomas-Fermi
out.a_TP_DFT = 2*exp(gammaE - 1+(2*Z*out.I1 + out.N*(1 - out.I2))/(Z+Z0))/alpha;%Tseng-Pratt
%
if strfind(model,'Thomas-Fermi')
    %
    a_DFT = out.a_TF_DFT;
    %
elseif strfind(model,'Tseng-Pratt')
    %
    a_DFT = out.a_TP_DFT;
    %
else
    a_DFT = NaN;
    out = NaN;
end

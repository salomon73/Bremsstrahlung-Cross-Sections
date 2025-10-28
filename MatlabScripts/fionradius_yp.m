function [Z0,N,q,a,nff,A,B,r0n,rn,bbar,dens,I1,I2,output] = fionradius_yp(Z,Z0,model, nY)
%
%   Calculation of the ion radius for any positive ion of an element of atomic number Z. 
%   It is calculated using  the density functional theory (DFT) by
%   different models. From the ion radius, it is possible to calculate the
%   screening function in Fokker-Planck equation as well as the screening
%   effect on the bremsstrahlung emission
%
%   Input:
%
%       - Z : Atomic number of the considered element 
%       - Z0 : charge of the ionized atom Z0 = Z - N
%       - model : DFT model ('Thomas-Fermi', 'Thomas-Fermi-Kirillov','Tseng-Pratt-Avdonina-Lamoureux','Tseng-Pratt-Botto','Tseng-Pratt-Thomas-Fermi','Full_screening',...)
%       - nY : number of Yukawa exponentials, for Multi-Yukawa only
%   Output:
%
%       - Z0
%       - N : Number of bound electrons (Z - Z0)
%       - q : ionization factor, q = 1 - N/Z
%       - a : output structure containing the the effective ion screening length (bohr radius units), model, and electron distribution
%       - nff : power law of the form factor
%       - A : Yukawa weight (only for multi-Yukawa)
%       - B : 
%       - r0n : 
%       - rn : 
%       - bbar :
%       - dens : 
%       - I1 :
%       - I2 :
%       - output : output structure of the numerical solution
%
%
%   Dependencies and toolboxes:
%       
%      Dependencies: 
%
%       - MYi_Zj.m                 % DFT multi-Yukawa coefficients
%       - load_yukawa.m            % Load multi-Yukawa data
%       - pc_dke_yp.m              % physical constants    
%       - simps.m                  % Simpson integration
%       - dft_integrals_yp.m       % Density functional theory            
%       - tf_potential_yp.m        % Thomas-Fermi potential model       
%       - trapz_dke_yp.m           % trapeze integration        
%           
%   Toolboxes:
%
%       - MATLAB
%
%
%   By Y. Peysson CEA-IRFM <yves.peysson@cea.fr>, A. Jardin IJF/PAN <axel.jardin@ifj.edu.pl>, Krzysztof Krol IFJ/PAN <krzysztof.krol@ifj.edu.pl>
%   and S. Guinchard <salomon.guinchard@epfl.ch>

if nargin < 2
    error('ERROR: not enough input parameters in fscreening_yp.m.');
end
%
if nargin == 2
    model = 'Thomas-Fermi';
    nY = 0; % dummy
end
%
[~,~,~,~,~,~,~,~,~,alpha,~,~,~,~,~,a0] = pc_dke_yp;
%
N = Z - Z0;%Number of bound electrons 
q = 1 - N/Z;%ionization factor
%
if strcmp(model,'No-screening') %fully stripped ion
    Z0 = Z;
    N  = 0;
    q  = 1;
    A = NaN;
end
%
if q < 1
    %
    if strcmp(model, 'Thomas-Fermi')
        %
        tstart = tic;
        [output.s0,B,output.err,output.s,output.y,output.exitflag] = tf_potential_yp(q);
        telapsed = toc(tstart);
        %
        output.telapsed = telapsed;%elapsed time
        %
        x0 = output.s0.^2;  %in b units where b = bbar*a0, bbar = (1/4)*(9*pi^2/(2*Z))^(1/3) and a0 is the Bohr radius
        x = output.s.^2;    %in b units where b = bbar*a0, bbar = (1/4)*(9*pi^2/(2*Z))^(1/3) and a0 is the Bohr radius
        %
        i0 = find(x==x0); %index where the potential vanishes
        %
        bbar = (1/4)*(9*pi^2/(2*Z))^(1/3);
        %
        C = Z/(4*pi)/bbar^3; %in a0 units
        dens = C*(output.y(1:i0,1)./x(1:i0)).^(3/2); %in b^(-3) units where a0 is the Bohr radius
        dens(1) = dens(2); %avoid Inf at x or rn = 0
        %
        if Z0 == Z
            r0n = 0; %in a0 units
        else
            r0n = bbar*x0; %in a0 units
        end
        rn = bbar*x(1:i0); %in a0 units
        rn(1) = eps;       %avoid NaN at x or rn = 0
        %
        [a,output.outint] = dft_integrals_yp(rn,dens,Z,N,Z0,model); %ion radius as calculated by L. Hesslow for fast electrons screening effect in magnetized plasmas
        %
        nff = 3/2; %for the form factor calculations
        %
        I1 = output.outint.I1;
        I2 = output.outint.I2;
        A = NaN;
        %
    elseif strcmp(model, 'Thomas-Fermi-Kirillov')
        %
        a = 2*(3*sqrt(pi))^(2/3)*N^(2/3)/Z/4/alpha;%in a0 units
        bbar = (1/4)*(9*pi^2/(2*Z))^(1/3);
        B = NaN;
        if Z0 == Z
            r0n = 0; %in a0 units
        else
            r0n = 10*((1/((1/Z)+q))^(1/3) - 1)*bbar; %in a0 units
        end
        rn = NaN; 
        dens = NaN;
        output = NaN;
        %
        nff = 3/2; %for the form factor calculations
        %
        I1 = NaN;
        I2 = NaN;
        A = NaN;
        %
    elseif strcmp(model, 'Multi-Yukawa')
        YukawaModel = load_yukawa('.',Z, Z0, nY, false);  % Load Yukawa data
        A = YukawaModel.A(Z0+1,:);                  % Extract Yukawa Weights 
        a = 2./alpha./YukawaModel.Lambda(Z0+1,:);   % Yukawa screening length          
        rn   =   NaN; 
        dens =   NaN;
        output = NaN;
        %
        nff = 2; %Form factor exponent
        %
        I1 = NaN;
        I2 = NaN;
    elseif strcmp(model, 'Tseng-Pratt-Avdonina-Lamoureux')
        %
        % N. B. Avdonina and R. H. Pratt, J. Quant. Spectrosc. Radiat. Transfer, 1993, 50, 4, pp. 349-358
        % M. Lamoureux and N. Avdonina, Physical Review E, 1997, 55, 1, pp. 912-926
        %
        nz = Z*(1/3 - 0.0020*Z);
        lambda0 = 0.8932*sqrt(Z);%a0^(-1) units with a0 Bohr radius  
        lambdaz2 = lambda0^2*(1.0 - (Z0/Z)^(nz+1))/(1.0 - (Z0/Z));
        lambdaz = sqrt(lambdaz2);
        %
        r = 0.99;
        n = 1000;
        h = (r - 1)/(r^(n-1) - 1); % 1st grid spacing
        x = [0, h*cumsum(r.^(0:(n-2)))];
        %
        if Z0 == Z
            rn = 0;
        else
            rn = 4*lambdaz*fliplr(1.0 - x)';%in a0 units
        end
        %
        dens = lambdaz2.*(Z-Z0).*exp(-lambdaz.*rn)./(4*pi*rn);
        test_dens = 4*pi*trapz_dke_yp([rn,dens.*rn.^2]);
        %
        B = NaN;
        r0n = NaN;
        bbar = NaN;
        %
        a = 2/lambdaz/alpha;
        %
        [output.a,output.outint] = dft_integrals_yp(rn,dens,Z,N,Z0,model);%ion radius as calculated by L. Hesslow for fast electrons screening effect in magnetized plasmas
        %
        nff = 2;%for the form factor calculations
        %
        I1 = output.outint.I1;
        I2 = output.outint.I2;
        A = NaN;
        %
    elseif strcmp(model,'Tseng-Pratt-Botto')
        %
        % D. J. Botto, J. McEnnan and R. H. Pratt, Phys. rev. A, Vol. 18, number 2 (1978) 580
        %
        nz = Z*(1/3 - 0.0020*Z);
        lambda0 = 0.9*Z^0.42;%a0^(-1) units with a0 Bohr radius  
        lambdaz2 = lambda0^2*(1.0 - (Z0/Z)^(nz+1))/(1.0 - (Z0/Z))^1.5;%the term .../(1.0 - (Z0/Z))^1.5 is added for an accurate dependencie with Z0/Z (notes Y. Peysson)
        %
        lambdaz = sqrt(lambdaz2);
        %
        r = 0.99;
        n = 1000;
        h = (r - 1)/(r^(n-1) - 1); % 1st grid spacing
        x = [0, h*cumsum(r.^(0:(n-2)))];
        %
        if Z0 == Z
            rn = 0;
        else
            rn = 4*lambdaz*fliplr(1.0 - x)';%in a0 units
        end
        %
        dens = lambdaz2.*(Z-Z0).*exp(-lambdaz.*rn)./(4*pi*rn);
        test_dens = 4*pi*trapz_dke_yp([rn,dens.*rn.^2]);
        %
        B = NaN;
        r0n = NaN;
        bbar = NaN;
        %
        a = 2/lambdaz/alpha;
        %
        [output.a,output.outint] = dft_integrals_yp(rn,dens,Z,N,Z0,model);%ion radius as calculated by L. Hesslow for fast electrons screening effect in magnetized plasmas
        %
        nff = 2;%for the form factor calculations
        %
        I1 = output.outint.I1;
        I2 = output.outint.I2;
        A = NaN;
        %
    elseif strcmp(model,'Tseng-Pratt-Thomas-Fermi')
        %
        % D. J. Botto, J. McEnnan and R. H. Pratt, Phys. rev. A, Vol. 18, number 2 (1978) 580
        %
        nz = Z*(1/3 - 0.0020*Z);
        lambda0 = 1.13*Z^(1/3);%a0^(-1) units with a0 Bohr radius (Thomas-Fermi value for neutral atom)
        lambdaz2 = lambda0^2*(1.0 - (Z0/Z)^(nz+1))/(1.0 - (Z0/Z))^1.5;%the term .../(1.0 - (Z0/Z))^1.5 is added for an accurate dependencie with Z0/Z (notes Y. Peysson)
        %
        lambdaz = sqrt(lambdaz2);
        %
        r = 0.99;
        n = 1000;
        h = (r - 1)/(r^(n-1) - 1); % 1st grid spacing
        x = [0, h*cumsum(r.^(0:(n-2)))];
        %
        if Z0 == Z
            rn = 0;
        else
            rn = 4*lambdaz*fliplr(1.0 - x)';%in a0 units
        end
        %
        dens = lambdaz2.*(Z-Z0).*exp(-lambdaz.*rn)./(4*pi*rn);
        test_dens = 4*pi*trapz_dke_yp([rn,dens.*rn.^2]);
        %
        B = NaN;
        r0n = NaN;
        bbar = NaN;
        %
        a = 2/lambdaz/alpha;
        %
        [output.a,output.outint] = dft_integrals_yp(rn,dens,Z,N,Z0,model);%ion radius as calculated by L. Hesslow for fast electrons screening effect in magnetized plasmas
        %
        nff = 2; %for the form factor calculations
        %
        I1 = output.outint.I1;
        I2 = output.outint.I2;
        A = NaN;
        %
    elseif strcmp(model,'Full-screening') || isempty(model)
        % 
        a = 0; %in a0 units
        bbar = NaN;
        B = NaN;
        r0n = NaN;
        rn = NaN; 
        dens = NaN;
        output = NaN;
        %
        nff = 2; %for the form factor calculations
        %
        I1 = NaN;
        I2 = NaN;
        A = NaN;
        %
    elseif strcmp(model,'DFT') || isempty(model)
        % 
        a = NaN; %in a0 units
        bbar = NaN;
        B = NaN;
        r0n = NaN;
        rn = NaN; 
        dens = NaN;
        output = NaN;
        %
        symbol = elements('atomic_number', Z, 'symbol');
        name = elements('atomic_number', Z, 'name');
        %
        load(['DFT/dft_',int2str(Z),'_',symbol,'_',name,'.mat'],'dft');
        %
        if isfield(dft,'I1')
            %
            I1 = dft.I1(Z0+1);
            I2 = dft.I2(Z0+1);
            %
        else
            error('Error: dft structure incomplete for fionradius_yp.m')
        end
        %
        nff = NaN; %for the form factor calculations
        %
        A = NaN;
    else
        error('ERROR: the atomic model is unknown.')
    end
    %
else
    %
    % No-screening model or fully stripped ion model
    %
    a = 0; %Ponctual charge
    bbar = NaN;
    B = NaN;
    r0n = NaN;
    rn = NaN; 
    dens = NaN;
    output = NaN;
    nff = 2;
    I1 = NaN;
    I2 = NaN;
    A = NaN;
    %
end
%
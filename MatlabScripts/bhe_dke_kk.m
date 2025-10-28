function [seBHE,seBH,Elwert,ec_out,k_out,c0_out,smodel,params] = bhe_dke_kk(ec_in,k_in,t0_in,Z,params,smodel)
%
%	Computes the e-i bremsstrahlung DDCS by integration of BH TDCS 
%   using an atomic form factor depending on the model. 
%	BH TDCS: see 1BS, 2BN formulas in H.W.Koch and J.W.Motz, Rev. Mod. Phys. 31, 4 (1959) 920
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
%		- ec_out: kinetic energy of the incoming electron (mc2) [p,m,n] 
%		- k_out: photon energy (mc2)  [p,m,n] 
%		- c0_out: cosine angle between the direction of displacement of the incoming electron and the photon emitted by bremsstrahlung (radian) [p,m,n]
%       - params = calculation parameters structure (for 1BS -> 2BN integration) 
%       - smodel = screening model structure
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
% by Y.Peysson CEA-DRFC <yves.peysson@cea.fr>, K.Krol PAN-IFJ <krzysztof.krol@ifj.edu.pl> and  A.Jardin PAN-IFJ <axel.jardin@ifj.edu.pl>
% and S. Guinchard <salomon.guinchard@epfl.ch>
% 
if nargin < 4
	infoyp(2,'Wrong number of input arguments for bhe_dke_kk');
	return;
end
%
if nargin < 5
    %
    params.nphi = 17;           % if NaN, 4-D calculations
    params.ntheta = 101;
    params.integration_method = 'trapz';
    %
	smodel.screening_model = ''; % no screening
    smodel.Z0 = '';              % fully screened ion charge
end
%
if nargin < 6
    %
    if isempty(params)
        params.nphi = 17;       % if NaN, 4-D calculations
        params.ntheta = 101;
        params.integration_method = 'trapz';
    end
    %
	smodel.screening_model = ''; % no screening
    smodel.Z0 = '';              % fully screened ion charge
end
%
[~,~,~,~,~,~,re,mc2,~,alpha] = pc_dke_yp; % Physics constant
%
% Scattered electron
%
if isfield(params,'nphi') && ~isnan(params.nphi) && ~isempty(params.nphi) 
    phi = linspace(0,2*pi,params.nphi);
end
%
theta  = relativistic_thetagrid_yp(max(ec_in), params.ntheta);
dtheta = theta(2) - theta(1);
c = cos(theta);
%
% For Elwert factor calculation (3-D) and some outputs
%
ec_out = repmat(ones(length(t0_in),1)*ec_in(:)'/mc2,[1,1,length(k_in)]);            % Relativistic units
k_out = shiftdim(repmat(ones(length(ec_in),1)*k_in(:)'/mc2,[1,1,length(t0_in)]),2); % Relativistic units
if (isscalar(ec_in)) && (isscalar(t0_in)) % scalar case
	k_out = reshape(k_out,1,1,length(k_in));
end
%
c0_out = repmat(t0_in(:)*ones(1,length(ec_in)),[1,1,length(k_in)]);
%
e0 = ec_out + 1;    %  Incoming electron (total) energy 
p0 = sqrt(e0.^2-1); %  memory saving
ep = e0 - k_out;    %  Initial electron momentum
p = sqrt(ep.^2-1);  %  final electron energy
%
mask = ec_out > k_out; % Photon are only emitted by electron of higher energies
%
ksi0 = alpha*Z*e0./p0;
ksi1 = alpha*Z*ep./p;
%
Elwert = (ksi1./ksi0).*((1-exp(-2*pi*ksi0))./(1-exp(-2*pi*ksi1))); % Elwert factor for Coulomb corrections
F0 = zeros(1,numel(smodel.Z0));
for iZ = 1:numel(smodel.Z0) % requested ionized states
    [smodel.Z0(iZ), ...
     smodel.N(iZ), ...
     smodel.qZ(iZ), ...
     smodel.a_Z0(iZ,:), ...
     smodel.nff(iZ), smodel.A_Z0(iZ,:)] = fionradius_yp(Z, smodel.Z0(iZ), smodel.screening_model, smodel.YKorder);

    F0(iZ) = Z - smodel.Z0(iZ);
end

%
if isfield(params,'nphi') && ~isnan(params.nphi) && ~isempty(params.nphi) 
    %
    % Extension to 5-D calculations for 1BS -> 2BN calculations
    %
    flag_mem = 0;
    %
    if strcmp(computer,'GLNXA64') || strcmp(computer,'MACI64') || strcmp(computer,'MACA64')
        %
        [~,memory] = system('sysctl hw.memsize | awk ''{print $2}''');
        memavail = str2double(strip(memory))/1e9; % in Gb
        disp(['Avail memory: ', num2str(memavail), 'Gb']);
        %
        allvars = whos;
        disp(['Memory used: ', num2str(sum([allvars.bytes])/1e9), 'Gb']); 
        %
        numblocks = 12; % # of blocks in the function
        memcall = numblocks*8*size(e0,1)*size(e0,2)*size(e0,3)*length(theta)*length(phi)/1e9; % in Gb 
        %
        if memcall > 2*memavail % takes into account of a possible swap memory (3xram)
            %
            disp(['WARNING: required memory (',num2str(memcall),' GBytes) is exceeding the maximum available RAM memory (',num22str(memavail),' GBytes).'])
            disp('Full block calculation is impossible on this computer with the required parameters.')
            disp('Loop in theta is performed.')
            %
            flag_mem = 1;
            %
        end
    else
        flag_mem = 1; % no memory check. Safe calculation mode considered by default
    end
    %
    if isfield(params,'flag_mem') && ~isempty(params.flag_mem) && (params.flag_mem ~= flag_mem)
        flag_mem = params.flag_mem; % overwrite automatic mode
        disp("override flag_mem")
    end
    %
    if flag_mem == 0 % full 5-D block calculations 
        %
        ec = repmat(repmat(ec_out,[1,1,1,length(theta)]),[1,1,1,1,length(phi)]);  % size(ec)
        k  = repmat(repmat(k_out,[1,1,1,length(theta)]),[1,1,1,1,length(phi)]);   % size(k)
        c0 = repmat(repmat(c0_out,[1,1,1,length(theta)]),[1,1,1,1,length(phi)]);  % size(c0)
        %
        c =   permute(repmat(ones(length(phi),1)*c(:)',[1,1,length(ec_in),length(k_in),length(t0_in)]),[5,3,4,2,1]);    %size(c)
        phi = permute(repmat(phi(:)*ones(1,length(theta)),[1,1,length(ec_in),length(k_in),length(t0_in)]),[5,3,4,2,1]); %size(phi)
        %
        s0 = sqrt(1-c0.^2); % sin(theta0)
        s = sqrt(1-c.^2);   % sin(theta)
        %
        e0 = ec + 1;        % Incoming electron (total) energy 
        clear ec;           % memory saving
        p0 = sqrt(e0.^2-1); % Initial electron momentum
        ep = e0 - k;        % final electron energy
        p = sqrt(ep.^2-1);  % Final electron momentum
        %
        q2 = p.^2 + p0.^2 + k.^2 - 2*p0.*k.*c0 + 2*p.*k.*c - 2*p0.*p.*(c.*c0 + s.*s0.*cos(phi)); % recoil momentum squared
        %
        b1 = p.^2.*s.^2./(ep - p.*c).^2.*(4*e0.^2 - q2);
        b0 = b1;
        clear b1;   %memory saving
        b2 = p0.^2.*s0.^2./(e0 - p0.*c0).^2.*(4*ep.^2 - q2);
        b0 = b0 + b2;
        clear b2;   %memory saving
        b3 = 2*p.*p0.*s.*s0.*cos(phi).*(4*ep.*e0-q2)./((ep - p.*c).*(e0 - p0.*c0));
        b0 = b0 - b3;
        clear b3;   %memory saving
        b4 = 2*k.^2.*(p.^2.*s.^2+p0.^2.*s0.^2-2*p.*p0.*s.*s0.*cos(phi))./((ep - p.*c).*(e0 - p0.*c0));
        b0 = b0 + b4;
        clear b4;   %memory saving
        %
        b0s = b0.*s;
        %
        if isfield(smodel,'Z0') && ~isempty(smodel.Z0) && isfield(smodel,'a_Z0') && ~isempty(smodel.a_Z0) && isfield(smodel,'nff') && ~isempty(smodel.nff)
            %
            ds3 = permute(repmat(zeros(length(smodel.Z0),1),[1,length(ec_in),length(k_in),length(t0_in)]),[4,2,3,1]); % allocate DDCS - size(ds3) 
            %
            for iZ = 1:length(smodel.Z0)
                % standardize shapes once 
                nIon = numel(smodel.Z0);
                if size(smodel.A_Z0,1) ~= nIon
                    smodel.A_Z0 = smodel.A_Z0.';    % use .' to avoid complex-conjugate transpose
                end
                if size(smodel.a_Z0,1) ~= nIon
                    smodel.a_Z0 = smodel.a_Z0.';
                end
                if ~isscalar(smodel.nff) && size(smodel.nff,1) ~= nIon
                    smodel.nff = smodel.nff.';
                end
                
                if strcmp(smodel.screening_model,'Multi-Yukawa')
                    sq = sqrt(q2);                 % same shape as q2
                    F = zeros(size(q2));
                
                    % take iZ row and make column vectors (nY x 1)
                    A_vec = smodel.A_Z0(iZ,:).';    % weights
                    a_vec = smodel.a_Z0(iZ,:).';    % ranges
                    if isscalar(smodel.nff)
                        nff_row = smodel.nff;       % scalar exponent
                    else
                        nff_row = smodel.nff(iZ,:); % 1 x nY
                    end
                
                    for kY = 1:numel(A_vec)
                        nff_k = nff_row(min(kY,end));   % handles scalar or per-term nff
                        F = F + (F0(iZ) * A_vec(kY)) ./ (1 + (sq .* (a_vec(kY)/2)).^nff_k);
                    end
                
                else
                    F = F0(iZ)./(1 + (sqrt(q2)*smodel.a_Z0(iZ)/2).^smodel.nff(iZ));
                end

                d = (Z - F).^2.*p./(k.*p0.*q2.^2);  % with screening effect described by the form factor (Born approximation)
                ds5 = d.*b0s;
                %
                if isfield(params,'integral') && strcmp(params.integral,'simps') && exist('simps', 'var')
                    ds4 = simps(linspace(0,2*pi,params.nphi),ds5,5);   %integration over phi
                    ds3(:,:,:,iZ) = simps(theta,ds4,4);                %integration over theta
                else
                    ds4 = trapz(linspace(0,2*pi,params.nphi),ds5,5);   %integration over phi
                    ds3(:,:,:,iZ) = trapz(theta,ds4,4);                %integration over theta
                end
            end
        else 
            % no screening
            d = Z^2.*p./(k.*p0.*q2.^2);% with screening effect described by the form factor (Born approximation)
            ds5 = d.*b0s;
            %
            if isfield(params,'integral') && strcmp(params.integral,'simps') && exist('simps', 'var')
                ds4 = simps(linspace(0,2*pi,params.nphi),ds5,5);  %integration over phi
                ds3 = simps(theta,ds4,4);                         %integration over theta
            else
                ds4 = trapz(linspace(0,2*pi,params.nphi),ds5,5);  %integration over phi
                ds3 = trapz(theta,ds4,4);                         %integration over theta
            end
            %
        end
        %
    else % 4-D block calculations (loop in theta)
        %
        if length(smodel.Z0) > 1
            ds3 = permute(repmat(zeros(length(smodel.Z0),1),[1,length(ec_in),length(k_in),length(t0_in)]),[4,2,3,1]);%size(ds3)
        end
        %
        ec = repmat(ec_out,[1,1,1,length(phi)]);    %size(ec)
        k = repmat(k_out,[1,1,1,length(phi)]);      %size(k)
        c0 = repmat(c0_out,[1,1,1,length(phi)]);    %size(c0)
        %
        s0 = sqrt(1-c0.^2);
        %
        cp = cos(phi);
        %
        cp = permute(repmat(cp(:),[1,length(ec_in),length(k_in),length(t0_in)]),[4,2,3,1]);%size(cp)
        %
        e0 = ec + 1;        %  Incoming electron (total) energy 
        clear ec;           %  memory saving
        p0 = sqrt(e0.^2-1); %  Initial electron momentum
        ep = e0 - k;        %  final electron energy
        p = sqrt(ep.^2-1);  %  Final electron momentum
        %
        for iZ = 1:length(smodel.Z0)
            %
            ds4 = 0;
            %
            for ic = 1:length(theta)
                %
                c = cos(theta(ic));
                s = sin(theta(ic));
                %
                q2 = p.^2 + p0.^2 + k.^2 - 2*p0.*k.*c0 + 2*p.*k.*c - 2*p0.*p.*(c.*c0 + s.*s0.*cp);
                %
                F = F0(iZ)./(1 + (sqrt(q2)*smodel.a_Z0(iZ)/2).^smodel.nff(iZ));%screening form factor
                %
                d = (Z - F).^2.*p./(k.*p0.*q2.^2); %with screening effect described by the form factor (Born approximation)
                %
                b1 = p.^2.*s.^2./(ep - p.*c).^2.*(4*e0.^2 - q2);
                b2 = p0.^2.*s0.^2./(e0 - p0.*c0).^2.*(4*ep.^2 - q2);
                b3 = 2*p.*p0.*s.*s0.*cp.*(4*ep.*e0-q2)./((ep - p.*c).*(e0 - p0.*c0));
                b4 = 2*k.^2.*(p.^2.*s.^2+p0.^2.*s0.^2-2*p.*p0.*s.*s0.*cp)./((ep - p.*c).*(e0 - p0.*c0));
                b0 = b1 + b2 - b3 + b4;
                %
                ds4 = ds4 + d.*b0.*s*dtheta;
                %
            end  
            %
            if isfield(params,'integral') && strcmp(params.integral,'simps') && exist('simps', 'var')
                if length(smodel.Z0) > 1
                    ds3(:,:,:,iZ) = simps(phi,ds4,4);  %integration over phi
                else 
                    ds3 = simps(phi,ds4,4);            %integration over phi
                end
            else
                if length(smodel.Z0) > 1
                    ds3(:,:,:,iZ) = trapz(phi,ds4,4);  %integration over phi
                else 
                    ds3 = trapz(phi,ds4,4);            %integration over phi
                end
            end
        end
        %
    end
    %
else
    %
    % Extension to 3-D calculations for 1BS -> 2BN calculations (integration in phi azimuthal angle has been performed analytically 4-D block calculations)
    %
    if length(smodel.Z0) > 1
            ds3 = permute(repmat(zeros(length(smodel.Z0),1),[1,length(ec_in),length(k_in),length(t0_in)]),[4,2,3,1]);%size(ds3)
    end
    %
    ec = repmat(ec_out,[1,1,1,length(theta)]);  %size(ec)
    k = repmat(k_out,[1,1,1,length(theta)]);    %size(k)
    c0 = repmat(c0_out,[1,1,1,length(theta)]);  %size(c0)
    %
    c = permute(repmat(c(:),[1,length(ec_in),length(k_in),length(t0_in)]),[4,2,3,1]);%size(c)
    %
    s0 = sqrt(1-c0.^2);
    s = sqrt(1-c.^2);
    %
    e0 = ec + 1;        %  Incoming electron (total) energy 
    p0 = sqrt(e0.^2-1); %  Initial electron momentum
    ep = e0 - k;        %  Final electron momentum
    p = sqrt(ep.^2-1);  %  final electron energy
    %
    aa =  -2*p0.*p.*s.*s0;
    bb = p.^2 + p0.^2 + k.^2 - 2*p0.*k.*c0 + 2*p.*k.*c - 2*p0.*p.*c.*c0;
    %
    d = Z.^2.*p./(k.*p0); % with screening effect described by the form factor (Born approximation)
    %
    b11 = p.^2.*s.^2./(ep - p.*c).^2;
    b21 = p0.^2.*s0.^2./(e0 - p0.*c0).^2;
    b31 = 2*p.*p0.*s.*s0./((ep - p.*c).*(e0 - p0.*c0));
    b41 = 2*k.^2./((ep - p.*c).*(e0 - p0.*c0));
    %
    for iZ = 1:length(smodel.Z0)
        %
        %ds4 = 0;
        %
        [J_1,J_2,J_3,J_4] = Jintegrals_tsengpratt(aa,bb,(smodel.a_Z0(iZ)/2).^2,smodel.qZ(iZ));
        %
        b0 = (4*e0.^2.*b11 + 4*ep.^2.*b21 + (p.^2.*s.^2 + p0.^2.*s0.^2).*b41).*J_1...
           - (b11 + b21).*J_2...
           - (4.*ep.*e0.*b31 + 2*p.*p0.*s.*s0.*b41).*J_3...
           + b31.*J_4;
        %
        %a1 = 4*e0.^2.*p.^2.*s.^2./(ep - p.*c).^2 + 4*ep.^2.*p0.^2.*s0.^2./(e0 - p0.*c0).^2 + 2*k.^2.*(p.^2.*s.^2 + p0.^2.*s0.^2)./((ep - p.*c).*(e0 - p0.*c0));
        %a2 = -(p.^2.*s.^2./(ep - p.*c).^2 + p0.^2.*s0.^2./(e0 - p0.*c0).^2);
        %a3 = -2*p.*p0.*s.*s0.*(4*ep.^2 + 4*e0.^2)./((ep - p.*c).*(e0 - p0.*c0));
        %a4 = 2*p.*p0.*s.*s0./((ep - p.*c).*(e0 - p0.*c0));
        %b0 = a1.*J_1 + a2.*J_2 + a3.*J_3 + a4.*J_4;
        %
        ds4 = d.*b0.*s;
        %
        if isfield(params,'integral') && strcmp(params.integral,'simps') && exist('simps', 'var')
            if length(smodel.Z0) > 1
                ds3(:,:,:,iZ) = simps(theta,ds4,4);%integration over theta
            else 
                ds3 = simps(theta,ds4,4);%integration over theta
            end
        else
            if length(smodel.Z0) > 1
                ds3(:,:,:,iZ) = trapz(theta,ds4,4);%integration over theta
            else 
                ds3 = trapz(theta,ds4,4);%integration over theta
            end
        end
    end
    %
end
%
if length(size(ds3)) > 3
    %
    mask = repmat(mask,[1,1,1,length(smodel.Z0)]);      %size(mask)
    Elwert = repmat(Elwert,[1,1,1,length(smodel.Z0)]);  %size(mask)
    %
    seBH = alpha*(re/(2*pi))^2.*ds3.*mask;
    %
    seBHE = seBH.*Elwert;
    %
else
    %
    seBH = alpha*(re/(2*pi))^2.*ds3.*mask;
    %
    seBHE = seBH.*Elwert;
    %    
end
%
seBHE(isnan(seBHE)) = 0;
%
seBHE(~isreal(seBHE)) = 0;
%
function [J_1,J_2,J_3,J_4] = Jintegrals_tsengpratt(a,b,c,qZ)
    %
    J_11 = zeros(size(a));
    J_21 = zeros(size(a));
    J_31 = zeros(size(a));
    J_41 = zeros(size(a));
    J_12 = zeros(size(a));
    J_13 = zeros(size(a));
    J_22 = zeros(size(a));
    J_23 = zeros(size(a));
    J_32 = zeros(size(a));
    J_33 = zeros(size(a));
    J_42 = zeros(size(a));
    J_43 = zeros(size(a));
    %
    delta = zeros(size(a));
    num = zeros(size(a));
    %
    A0 = zeros(size(a));
    B0 = zeros(size(a));
    C0 = zeros(size(a));
    D0 = zeros(size(a));
    A = zeros(size(a));
    B = zeros(size(a));
    T1 = zeros(size(a));
    T2 = zeros(size(a));
    T3 = zeros(size(a));
    T4 = zeros(size(a));
    %
    i0 = find(a == 0);
    in0 = find(a ~= 0);
    %
    delta(in0) = b(in0)./a(in0);% delta < -1
    num(in0) = sqrt(delta(in0).^2 - 1);
    %
    J_11(in0) = -2*pi*delta(in0)./num(in0).^3./a(in0).^2;
    %
    J_21(in0) = -2*pi./num(in0)./a(in0);
    %
    J_31(in0) = 2*pi./num(in0).^3./a(in0).^2;
    %  
    J_41(in0) = 2*pi*(1 + delta(in0)./num(in0))./a(in0);
    %
    J_11(i0) = 2*pi./b(i0).^2;
    J_21(i0) = 2*pi./b(i0);
    J_31(i0) = J_31(i0)*0;
    J_41(i0) = J_41(i0)*0;
    %
    if isnan(c) && isnan(qZ) %no screening
        %
    else %Tseng Pratt screening model
        %
        A0(in0) = delta(in0) - 1;               % A0 < 0
        B0(in0) = delta(in0) + 1;               % B0 < 0
        %
        C0(in0) = 1 + a(in0).*c.*(delta(in0) - 1);   % C0 > 0
        D0(in0) = 1 + a(in0).*c.*(delta(in0) + 1);   % D0 > 0
        %
        A(in0) = B0(in0)./A0(in0); % A > 0
        B(in0) = D0(in0)./C0(in0); % B > 0
        %
        T1(in0) = (A(in0).^2-1+2*B(in0).*(1-A(in0)))./(B(in0)-A(in0)).^2;
        T2(in0) = (1-A(in0)).^2./(B(in0)-A(in0));
        T3(in0) = (B(in0)-1).^2./(B(in0)-A(in0)).^2;
        %
        J_12(in0) = 2*pi./a(in0).^2.*1./(A0(in0).^2.*C0(in0)).*(T1(in0)./sqrt(A(in0))+T2(in0)./(2*A(in0).^(3/2))+T3(in0)./sqrt(B(in0)));
        %
        T1(in0) = (3*A(in0).^2.*B(in0)+3*B(in0)+3*A(in0)-6*A(in0).*B(in0)-A(in0).^3-2)./(B(in0)-A(in0)).^3;
        T2(in0) = (1-A(in0)).^3./(B(in0)-A(in0)).^2;
        T3(in0) = (B(in0).^3-3*A(in0).*B(in0).^2-3*B(in0)-3*A(in0)+6*A(in0).*B(in0)+2)./(B(in0)-A(in0)).^3;
        T4(in0) = (3*B(in0).^2-B(in0).^3-3*B(in0)+1) ./ (B(in0)-A(in0)).^2;
        %
        J_13(in0) = 2*pi./a(in0).^2.*1./(A0(in0).^2.*C0(in0).^2).*(T1(in0)./sqrt(A(in0))+T2(in0)./(2*A(in0).^(3/2))+T3(in0)./sqrt(B(in0))+T4(in0)./(2*B(in0).^(3/2)));
        %
        J_22(in0) = 2*pi./a(in0).*1./(A0(in0).*C0(in0)).*1./(B(in0)-A(in0)).*((1-A(in0))./sqrt(A(in0))+(B(in0)-1)./sqrt(B(in0)));
        %
        T1(in0) = (B(in0).^2-1+2*A(in0).*(1-B(in0)))./(A(in0)-B(in0)).^2;
        T2(in0) = (1-B(in0)).^2./(A(in0)-B(in0));
        T3(in0) = (A(in0)-1).^2./(A(in0)-B(in0)).^2;
        %
        J_23(in0) = 2*pi./a(in0).*1./(A0(in0).*C0(in0).^2).*(T1(in0)./sqrt(B(in0))+T2(in0)./(2*B(in0).^(3/2))+T3(in0)./sqrt(A(in0)));
        %
        J_32(in0) = 1./a(in0).*J_22(in0)-delta(in0).*J_12(in0);
        %
        J_33(in0) = 1./a(in0).*J_23(in0)-delta(in0).*J_13(in0);
        %
        J_42(in0) = 2*pi./a(in0).*1./sqrt(C0(in0).*D0(in0))-delta(in0).*J_22(in0);
        %
        J_43(in0) = pi./a(in0).*(B(in0)+1)./(B(in0).^(3/2).*C0(in0).^2)-delta(in0).* J_23(in0);
        %
        J_12(i0) = 2*pi./b(i0).^2./(1 + b(i0).*c);
        J_13(i0) = 2*pi./b(i0).^2./(1 + b(i0).*c).^2;
        J_22(i0) = 2*pi./b(i0)./(1 + b(i0).*c);
        J_23(i0) = 2*pi./b(i0)./(1 + b(i0).*c).^2;
        J_32(i0) = 0*J_32(i0);
        J_33(i0) = 0*J_33(i0);
        J_42(i0) = 0*J_42(i0);
        J_43(i0) = 0*J_43(i0);
        %
    end
    % 
    if isnan(qZ)
        %
        J_1 = J_11;
        J_2 = J_21;
        J_3 = J_31;
        J_4 = J_41;
        %
    else
        %
        J_1 = J_11 - 2*(1-qZ).*J_12 + (1-qZ).^2.*J_13;
        J_2 = J_21 - 2*(1-qZ).*J_22 + (1-qZ).^2.*J_23;
        J_3 = J_31 - 2*(1-qZ).*J_32 + (1-qZ).^2.*J_33;
        J_4 = J_41 - 2*(1-qZ).*J_42 + (1-qZ).^2.*J_43;
        %
    end
    %
    %  [err,output] = test_J_screening_integrals(squeeze(a(1,1,1,1:10)),squeeze(b(1,1,1,1:10)),squeeze(c(1,1,1,1)));%for test only
    %
end
%





end

function [s0,B,err,s,y,exitflag] = tf_potential_yp(q,smax)
    %
    %   Calculation of the solution of the dimensionless Thomas-Fermi equation for any screening factor q, by solving two coupled ordinary differential equations. 
    %
    %   INPUTS:
    %
    %       - q : screening factor (1 - N/Z) [1,1]
    %       - smax : upper limit of the radial grid (normalized to Thomas-Fermi length) [1,1]
    %
    %   OUTPUTS:
    %
    %       - s = normalized modified radial coordinate up to value where dimensionless potential vanishes (positive ion only), x = sqrt(s). 
    %       - y = [w(s) = phi(x),u(s) = w'(s)/s] dimensionless Thomas-Fermi equation
    %       - B = derivative at the origin
    %       - err = qnum - qref for determining b_qref consistently
    %       - exitflag = zero search output flag
    %
    %   By Y. Peysson CEA-IRFM <yves.peysson@cea.fr>, A. Jardin IJF/PAN <axel.jardin@ifj.edu.pl>, Krzysztof Krol IFJ/PAN <krzysztof.krol@ifj.edu.pl>
    %
    if nargin == 0 
        %
        q = 0.5;%weakly ionized atom
        %
    end
    %
    [B_init,x0] = tf_init(q);
    %
    if nargin < 2
        %
        smax = sqrt(x0)*1.05;%5% margin in the search domain
        %
    end
    %
    fun_tf = @(B) tf_root(B,q,smax);
    %
    [B,~,exitflag,output] = fzero(fun_tf,B_init);%Find the first order derivative B for a given q value by solving the dimensionless Thomas-Fermi equation
    %
    [err,y,s,s0] = tf_root(B,q,smax);%Calculate the dimensionless potential with B and q
    %
    function [err,y,s,s0] = tf_root(B,q,smax)
        %
        sspan = [0,smax];
        %
        y1_0 = 1;                   %w@s=0
        y2_0 = -2*B;%u@s=0
        %
        options = odeset('MaxStep',5e-4*smax^2);%adaptative step size for accurate calculation of the number of bounded electrons
        %
        [s,y] = ode45( @tf_equations ,sspan ,[y1_0,y2_0],options);%evaluates r.h.s. of the two coupled non-linear odes
        %
        y = real(y);%remove useless imaginary part
        %
        imax = find(y(:,1)>=0, 1,'last');%index of the value that is closest to zero (vanishing potential)
        err = -y(imax,2)*s(imax)^2/2 - q;
        s0 = s(imax);
        %
    end
    %
    function [B,x0] = tf_init(q)
        %
        % From "Semiclassical Theory of Atoms" by Berthold-Georg Englert (Chapter II),  Lecture Notes in Physics, Springer-Verlag, Germany, 1988
        % 
        nsz_ref = [0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95];
        x0_ref = [0.416269,0.685790,0.934348,1.179253,1.428919,1.689292,1.965691,2.263681,2.589715,2.951825,3.360561,3.830452,4.382486,5.048683,5.881272,6.973385,8.513784,10.92728,16.10273];
        B_ref = [3.020996,2.233243,1.952470,1.813524,1.734116,1.684993,1.653119,1.631819,1.617337,1.607410,1.600602,1.595965,1.592853,1.590815,1.589530,1.588763,1.588345,1.588149,1.588081];
        %
        x0 = interp1(nsz_ref,x0_ref,1-q,'spline','extrap');
        B = interp1(nsz_ref,B_ref,1-q,'spline','extrap');
        %
    end
    %
    function dyds = tf_equations(s,y) % Dimensionless Thomas-Fermi equation transformed in two coupled non-linear odes
        %
        dyds = zeros(2,1);%(w,u)
        %
        dyds(1) = s*y(2);%w'
        dyds(2) = 4*y(1)^1.5;%u'
        %
     end
     %
end
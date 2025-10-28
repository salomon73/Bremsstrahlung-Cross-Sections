function [stat_error_abs, syst_error_abs] = errors_cross_sec(E0,xp,yp, Motz)
%
% Systematic and statistical errors on the data given in input xp, yp,
% corresponding to emitted photon energy (xp) and measured differential
% cross section (yp). 
%
% Refs : 
% W. Koch and J. W. Motz, Bremsstrahlung Cross-Section Formulas and Related Data, Rev. Mod. Phys. , 1959, 31, 4, pp. 920
% N. Starfelt and H. W. Koch, Differential Cross-Section Measurements of Thin-Target Bremsstrahlung Produced by 2.7- to 9.7-Mev Electrons, Phys. Rev., 1956, 102, 6, pp. 1598-1612
% D. H. Rester and W. E. Dance, Bremsstrahlung Cross-Section Measurements at Incident Electron Energies of 1.0, 1.7, and 2.5 MeV
% J. W. Motz, Bremsstrahlung Differential Cross-Section Measurements for 0.5- and 1.0-Mev Electrons
%
%   Input:
%
%       - E0: incoming electron energy (MeV)
%       - xp: emitted photon energies (MeV)
%       - yp: differential cross-sec (cm^2/MeV/sr)
%
%	Output: 
%
%       - stat_error_abs: statistical errors (quadratically added)
%       - syst_error_abs: systematic errors (linearly added)
%
%   WARNING : energy of incoming electrons are in MeV. 
%
%   S. Guinchard EPFL <salomon.guinchard@epfl.ch>
%

if nargin < 4
    Motz = false;
end

    k_over_E0 = xp / E0;
    k_over_E0_table = [0.1, 0.4, 0.8, 0.95];
    
    % errors from Starfelt - Koch
    if E0 == 9.66
        % Table II from Starfelt and Koch for 9.66 MeV
        stat_error_percent_table = [3, 5, 6.5, 17.5]; 
        
        % Interpolate statistical error 
        stat_error_percent = interp1(k_over_E0_table, stat_error_percent_table, k_over_E0, 'linear', 'extrap');
        
        % Absolute statistical error
        stat_error_abs = yp .* (stat_error_percent / 100);
        
        % Systematic error 
        syst_error_percent = 3+4+1.5;  % see group b first column in Table II.
        syst_error_abs = yp .* (syst_error_percent / 100);
    elseif E0 == 4.54
        % Table II from Starfelt and Koch for 4.54 MeV
        stat_error_percent_table = [4, 6, 6.5, 14.5]; 
        
        % Interpolate statistical error 
        stat_error_percent = interp1(k_over_E0_table, stat_error_percent_table, k_over_E0, 'linear', 'extrap');
        
        % Absolute statistical error
        stat_error_abs = yp .* (stat_error_percent / 100);
        
        % Systematic error 
        syst_error_percent = 6+5;  % see group b second column in Table II.
        syst_error_abs = yp .* (syst_error_percent / 100);
    elseif E0 == 2.72
        % Table II from Starfelt and Koch for 2.72 MeV
        stat_error_percent_table = [6, 8.5, 19.0, 26.0]; 
        
        % Interpolate statistical error 
        stat_error_percent = interp1(k_over_E0_table, stat_error_percent_table, k_over_E0, 'linear', 'extrap');
        
        % Absolute statistical error
        stat_error_abs = yp .* (stat_error_percent / 100);
        
        % Systematic error 
        syst_error_percent = 2;  % see group b third column in Table II.
        syst_error_abs = yp .* (syst_error_percent / 100);

    % Errors from Rester - Dance
    elseif E0 == 2.50
        % systematic (linearly added): ±2% (current) + ±5% (thickness)
        syst_error_percent = 2 + 5;
        syst_error_abs = yp .* (syst_error_percent / 100);
    
        % statistical & angle & response (added in quadrature)
        % Angle: 12% 
        angle_err = 12;
    
        % Statistical: < 8% below 0.5 E0, up to 30% at 0.9 E0
        stat_err = interp1([0.5, 0.9], [8, 30], k_over_E0, 'linear', 'extrap');
    
        % Spectrometer: < 3% below 0.9 E0, 30–50% above
        response_err = zeros(size(k_over_E0));
        response_err(k_over_E0 < 0.9) = 3;
        response_err(k_over_E0 >= 0.9) = interp1([0.9, 1.0], [30, 50], k_over_E0(k_over_E0 >= 0.9), 'linear', 'extrap');
    
        % Total point-by-point error (RMS)
        stat_error_percent = sqrt(angle_err.^2 + stat_err.^2 + response_err.^2);
        stat_error_abs = yp .* (stat_error_percent / 100);
    
    elseif E0 == 1.70
        syst_error_percent = 2 + 5;
        syst_error_abs = yp .* (syst_error_percent / 100);
    
        angle_err = 9;  % interpolated approx between 6% and 12%
        stat_err = interp1([0.5, 0.9], [6, 25], k_over_E0, 'linear', 'extrap');
        response_err = zeros(size(k_over_E0));
        response_err(k_over_E0 < 0.9) = 3;
        response_err(k_over_E0 >= 0.9) = interp1([0.9, 1.0], [30, 50], k_over_E0(k_over_E0 >= 0.9), 'linear', 'extrap');
    
        stat_error_percent = sqrt(angle_err.^2 + stat_err.^2 + response_err.^2);
        stat_error_abs = yp .* (stat_error_percent / 100);
    
    elseif E0 == 1.00
        syst_error_percent = 2 + 5;
        syst_error_abs = yp .* (syst_error_percent / 100);
    
        angle_err = 6;
        stat_err = interp1([0.5, 0.9], [4, 15], k_over_E0, 'linear', 'extrap');
        response_err = zeros(size(k_over_E0));
        response_err(k_over_E0 < 0.9) = 3;
        response_err(k_over_E0 >= 0.9) = interp1([0.9, 1.0], [30, 50], k_over_E0(k_over_E0 >= 0.9), 'linear', 'extrap');
    
        stat_error_percent = sqrt(angle_err.^2 + stat_err.^2 + response_err.^2);
        stat_error_abs = yp .* (stat_error_percent / 100);
    
    elseif E0==0.2 
        % Rester & Edmonson (1972) at 0.2 MeV
        % Average error:
        % - 5% below 100 keV
        % - 10% from 100 to 180 keV
        % - 30% above 180 keV
        
        syst_error_abs = zeros(size(yp));  % not separately identified
        stat_error_abs = zeros(size(yp));
        
        for i = 1:length(xp)
            if xp(i) < 0.1
                total_error_percent = 5;
            elseif xp(i) >= 0.1 && xp(i) < 0.18
                total_error_percent = 10;
            else
                total_error_percent = 30;
            end
            stat_error_abs(i) = yp(i) * total_error_percent / 100;
        end
    elseif E0==0.18 

        %  Al 180 keV Klasmeier data
        syst_error_abs = zeros(size(yp));  
        stat_error_abs = zeros(size(yp));

    
    % Errors from Motz
    elseif (E0 == 1.00 && Motz)
        % systematic
        % Photon angle + current + uniformity + beam-target angle + flux
        syst_error_percent = 10 + 2 + 3 + 2 + 5;
        syst_error_abs = yp .* (syst_error_percent / 100);
    
        % statistical
        stat_err = zeros(size(k_over_E0));
        stat_err(k_over_E0 < 0.8) = 5;
        stat_err(k_over_E0 >= 0.8 & k_over_E0 < 0.9) = 20;
        stat_err(k_over_E0 >= 0.9) = 30;
    
        % Add in quadrature:
        window_err = 4;
        graph_integ_err = 3;
    
        stat_error_percent = sqrt(stat_err.^2 + window_err^2 + graph_integ_err^2);
        stat_error_abs = yp .* (stat_error_percent / 100);

    elseif E0 == 0.50 
        % systematic
        syst_error_percent = 10 + 2 + 3 + 2 + 5;  % From Motz summary
        syst_error_abs = yp .* (syst_error_percent / 100);
    
        % statistical
        stat_err = zeros(size(k_over_E0));
        stat_err(k_over_E0 < 0.8) = 5;
        stat_err(k_over_E0 >= 0.8 & k_over_E0 < 0.9) = 20;
        stat_err(k_over_E0 >= 0.9) = 30;
    
        window_err = 4;
        graph_integ_err = 3;
    
        stat_error_percent = sqrt(stat_err.^2 + window_err^2 + graph_integ_err^2);
        stat_error_abs = yp .* (stat_error_percent / 100);

    else 
         %error('ERROR: energy should be in [0.50, 1.00, 1.70, 2.50, 2.72, 4.54, 9.66] MeV.')
         syst_error_abs = zeros(1,numel(xp));
         stat_error_abs = zeros(1,numel(xp));
    end
end
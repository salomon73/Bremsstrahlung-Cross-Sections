function v = kinetic_energy_to_velocity(Ek_keV)
    % Constants
    c = 299792458;               % speed of light [m/s]
    mc2 = 511;                   % electron rest energy [keV]
    
    gamma = Ek_keV / mc2 + 1;
    v = c * sqrt(1 - 1 ./ gamma.^2);  % [m/s]
end


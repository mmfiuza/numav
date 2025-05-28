module Analytic

export simulate_analytic_rectangular_room

function calculate_modes(position, lx, ly, lz, speed_of_sound, max_mode_idx)

    x = position[1];
    y = position[2];
    z = position[3];
    psi = Vector{Float64}(undef, (max_mode_idx+1)*(max_mode_idx+1)*(max_mode_idx+1));
    freq_modes = Vector{Float64}(undef, (max_mode_idx+1)*(max_mode_idx+1)*(max_mode_idx+1));
    i = 0;
    for nz in 0:max_mode_idx
        for ny in 0:max_mode_idx
            for nx in 0:max_mode_idx
                i = i+1;
                freq_modes[i] = speed_of_sound/2 * sqrt((nx/lx)^2 + (ny/ly)^2 + (nz/lz)^2);
                psi[i] = cos(nx*pi*x/lx) .* cos(ny*pi*y/ly) .* cos(nz*pi*z/lz);
            end
        end
    end
    sort_order = sortperm(freq_modes);
    freq_modes = freq_modes[sort_order];
    psi = psi[sort_order];

    return psi, freq_modes
end

function simulate_analytic_rectangular_room(room_disensions, source_position, receiver_position, 
temperature_celsius, atmospheric_pressure, freq_vector, volume_velocity_amplitude::Function)

    rho_air = atmospheric_pressure / (287*(temperature_celsius + 273.15)); # [kg/m^3]
    c_air = sqrt(1.400 * 287.0 * (temperature_celsius + 273.15)); # [m/s]

    lx = room_disensions[1]
    ly = room_disensions[2]
    lz = room_disensions[3]
    volume = lx*ly*lz;

    # calculate modes
    max_mode_idx = 100
    psi_source, freq_modes = calculate_modes(source_position, lx, ly, lz, c_air, max_mode_idx)
    psi_receiver, _ = calculate_modes(receiver_position, lx, ly, lz, c_air, max_mode_idx)
    omega_modes = 2pi*freq_modes
    k_modes = omega_modes/c_air

    omega = 2pi*freq_vector
    k = omega/c_air

    complex_pressure_amplitude = 
        im*omega*rho_air.*volume_velocity_amplitude.(freq_vector)/volume .* sum(
        transpose(psi_receiver).*transpose(psi_source) ./ (k.^2 .- transpose(k_modes).^2),
        dims=2
    )
    return complex_pressure_amplitude
end

end
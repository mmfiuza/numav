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
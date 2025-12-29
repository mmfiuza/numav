# Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>
#
# READ_NMVR Read NUMAV Result (.nmvr) file format
#   data = read_nmvr(filename) reads the .nmvr file and returns a
#   structure with all the organized data.

function read_nmvr(filename::String)

    # check file extension
    if !endswith(filename, ".nmvr")
        @warn "The file name does not end with .nmvr"
    end

    # open file for reading (binary, little-endian)
    fid = open(filename, "r")
    try
        # check signature
        signature = String(read(fid, 8))
        if signature != "numavrst"
            error("The file selected is not written in the nmvr format.")
        end

        # read sim_type chunk
        chunk_id = String(read(fid, 8))
        if chunk_id != "sim_type"
            error("Simulation type not found.")
        end
        
        chunk_size = read(fid, UInt64)
        if chunk_size == 40
            sim_type = strip(String(read(fid, 32)))
            if sim_type == "acousticfem_____freq_dom3_dimens"
                return ac_fem_freq_d3(fid)
            else
                error("Simulation type not found.")
            end
        else
            error("Simulation type not found.")
        end
    finally
        close(fid)
    end
end

function ac_fem_freq_d3(fid::IOStream)

    # Initialize structure
    data = Dict{String, Any}()
    data["simulation_type"] = Dict{String, String}(
        "phenomenon" => "acoustic",
        "numerical_method" => "fem_____",
        "domain" => "freq_dom",
        "dimension" => "3_dimens"
    )
    data["freq_stp"] = Float64[]
    data["ni2coord"] = Matrix{Float64}(undef, 0, 0)
    data["sei_2_ni"] = Matrix{UInt64}(undef, 0, 0)
    data["vei_2_ni"] = Matrix{UInt64}(undef, 0, 0)
    data["cpx_pres"] = ComplexF64[]
    
    # read element order
    order = strip(String(read(fid, 8)))
    if contains(order, "1st_ord_")
        element_order = 1
        data["simulation_type"]["element_order"] = "1st_ord_"
    elseif contains(order, "2nd_ord_")
        element_order = 2
        data["simulation_type"]["element_order"] = "2nd_ord_"
    else
        error("Unknown element order: $order")
    end
    
    # Read chunks until end of file
    while !eof(fid)
        
        chunk_id = String(read(fid, 8))
        if chunk_id == "numaveof"
            break
        end
        chunk_size = read(fid, UInt64)
    
        if chunk_id == "freq_stp"
            n_elements = div(chunk_size, 8)
            data["freq_stp"] = [read(fid, Float64) for _ in 1:n_elements]
            
        elseif chunk_id == "ni2coord"
            coordinate_count = div(chunk_size, 8)
            if mod(coordinate_count, 3) != 0
                error("ni2coord data length is not a multiple of 3")
            end
            coords = [read(fid, Float64) for _ in 1:coordinate_count]
            node_count = div(coordinate_count, 3)
            data["ni2coord"] = reshape(coords, (3, node_count))
            
        elseif chunk_id == "sei_2_ni"
            ni_count = div(chunk_size, 8)
            # reshape based on element order
            if element_order == 1
                nodes_per_element = 3
            else
                nodes_per_element = 6
            end
            if mod(ni_count, nodes_per_element) != 0
                error("sei_2_ni data length is not a multiple of $nodes_per_element.")
            end
            
            indices = [read(fid, UInt64) for _ in 1:ni_count]
            indices = Int.(indices)
            element_count = div(ni_count, nodes_per_element)
            data["sei_2_ni"] = 
                reshape(indices, (nodes_per_element, element_count))
            
        elseif chunk_id == "vei_2_ni"
            ni_count = div(chunk_size, 8)
            # reshape based on element order
            if element_order == 1
                nodes_per_element = 4
            else
                nodes_per_element = 10
            end
            if mod(ni_count, nodes_per_element) != 0
                error("vei_2_ni data length is not a multiple of $nodes_per_element.")
            end
            
            indices = [read(fid, UInt64) for _ in 1:ni_count]
            indices = Int.(indices)
            element_count = div(ni_count, nodes_per_element)
            data["vei_2_ni"] =
                reshape(indices, (nodes_per_element, element_count))
        
        elseif chunk_id == "cpx_pres"
            complex_count = div(chunk_size, 16)
            
            # read real parts (skip 8 bytes between reads)
            real_parts = Vector{Float64}(undef, complex_count)
            for i in 1:complex_count
                real_parts[i] = read(fid, Float64)
                skip(fid, 8)
            end
            
            # go back to beginning of complex data
            seek(fid, position(fid) - chunk_size)
            
            # read imaginary parts (skip 8 bytes between reads)
            imag_parts = Vector{Float64}(undef, complex_count)
            skip(fid, 8)
            for i in 1:complex_count
                imag_parts[i] = read(fid, Float64)
                skip(fid, 8)
            end
            
            # adjust position for next read
            skip(fid, -8)
            
            data["cpx_pres"] = complex.(real_parts, imag_parts)

        else
            # unknown chunk type - skip it
            @warn "Unknown chunk ID: $chunk_id. Skipping $chunk_size bytes."
            skip(fid, chunk_size)
        end
    end
    
    # validate data consistency
    if isempty(data["freq_stp"])
        error("Could not determine frequency steps.")
    end
    if isempty(data["ni2coord"])
        error("Could not determine nodes.")
    end
    if isempty(data["sei_2_ni"])
        error("Could not determine surface elements.")
    end
    if isempty(data["vei_2_ni"])
        error("Could not determine volume elements.")
    end
    if isempty(data["cpx_pres"])
        error("Could not determine pressure solution.")
    end

    # reshape complex pressure to (ni_count × freq_count)
    node_count = size(data["ni2coord"], 2)
    freq_count = length(data["freq_stp"])
    complex_count = length(data["cpx_pres"])
    if node_count * freq_count != complex_count
        error("cpx_pres dimensions mismatch. Expected $node_count×$freq_count=$(node_count*freq_count), got $complex_count elements.")
    end
    data["cpx_pres"] = reshape(data["cpx_pres"], (node_count, freq_count))

    # convert to one-based index
    data["sei_2_ni"] = data["sei_2_ni"] .+ 1
    data["vei_2_ni"] = data["vei_2_ni"] .+ 1

    return data
end
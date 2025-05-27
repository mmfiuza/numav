using Printf

mutable struct Mesh
    nodes::Matrix{Float64}
    elem_surf::Matrix{Int}
    elem_vol::Matrix{Int}
    id_surf::Vector{Int}
    id_vol::Vector{Int}
    num_nodes::Int
    num_elem_surf::Int
    num_elem_vol::Int
end

function read_bdf(infilename)

    function replace_pm(S_in)
        S = S_in
        positions = vcat(findall(==('+'), S), findall(==('-'), S))
        sort!(positions, rev=true)
        
        for pos in positions
            if pos != 1 && S[prevind(S, pos)] in '0':'9'
                S = S[1:prevind(S, pos)] * "e" * S[pos:end]
            end
        end
        return S
    end
    
    lines = readlines(infilename)
    
    # First pass: Counting
    num_nodes = 0
    num_ele_surf = 0
    num_ele_vol = 0
    for line in lines
        if(startswith(line, "GRID")) num_nodes+=1 end
        if(startswith(line, "CTRIA3")) num_ele_surf+=1 end
        if(startswith(line, "CTETRA")) num_ele_vol+=1 end
    end
    
    # Initialize arrays
    nodes = Matrix{Float64}(undef, 3, num_nodes)
    elements_surf = Matrix{Int}(undef, 3, num_ele_surf)
    elements_vol = Matrix{Int}(undef, 4, num_ele_vol)
    id_surf = Vector{Int}(undef, num_ele_surf)
    id_vol = Vector{Int}(undef, num_ele_vol)
    
    # Second pass: Parse data
    idx_cur_node = 0
    idx_cur_ele_surf = 0
    idx_cur_ele_vol = 0

    for line in lines
        if startswith(line, "GRID")
            idx_cur_node += 1
            if startswith(line, "GRID*")
                # Long format parsing
                onn = parse(Int, line[9:24])
                sx = replace_pm(line[41:56])
                sy = replace_pm(line[57:72])
                sz = replace_pm(line[73:88])
            else
                # Short format parsing
                onn = parse(Int, line[9:16])
                sx = replace_pm(line[25:32])
                sy = replace_pm(line[33:40])
                sz = replace_pm(line[41:48])
            end
            
            nodes[:,idx_cur_node] .= [parse(Float64,sx), parse(Float64,sy), parse(Float64,sz)]
        
        elseif startswith(line, "CTRIA3")
            idx_cur_ele_surf += 1
            data = parse.(Int, split(line[17:end]))
            id_surf[idx_cur_ele_surf] = data[1]
            elements_surf[:,idx_cur_ele_surf] = [data[2], data[3], data[4]]
        
        elseif startswith(line, "CTETRA")
            idx_cur_ele_vol += 1
            data = parse.(Int, split(line[17:end]))
            id_vol[idx_cur_ele_vol] = data[1]
            elements_vol[:,idx_cur_ele_vol] = [data[2], data[3], data[4], data[5]]
        
        end
    end
    
    return Mesh(
        nodes,
        elements_surf,
        elements_vol,
        id_surf,
        id_vol,
        num_nodes,
        num_ele_surf,
        num_ele_vol
    )
end
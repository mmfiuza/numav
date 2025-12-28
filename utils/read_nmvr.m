% Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>
%
% READ_NMVR Read NUMAV Result (.nmvr) file format
%   data = READ_NMVR(filename) reads the .nmvr file and returns a
%   structure with all the organized data.

function data = read_nmvr(filename)

last_5_chars = filename(end-4:end);
if ~strcmp(last_5_chars, '.nmvr')
    warning('The file name does not end with .nmvr')
end

% open file for reading (binary, little-endian)
fid = fopen(filename, 'r', 'l');
if fid == -1
    error('Cannot open file: %s', filename);
end

% check signature
signature = fread(fid, 8, 'char=>char')';
if ~strcmp(signature, 'numavrst')
    error('The file selected is not witten in the nmvr format.');
end

% read sim_type chunk
chunk_id = fread(fid, 8, 'char=>char')';
if ~strcmp(chunk_id, 'sim_type')
    error('Simulation type not found.')
end
chunk_size = fread(fid, 1, 'uint64=>uint64');
if chunk_size == 40
    sim_type = strtrim(fread(fid, 32, 'char=>char')');
    if strcmp(sim_type, 'acousticfem_____freq_dom3_dimens')
        ac_fem_freq_d3();
    else
        error('Simulation type not found.')
    end
else
    error('Simulation type not found.')
end

% close file
fclose(fid);

function ac_fem_freq_d3()

    % Initialize structure
    data.simulation_type.phenomenon = 'acoustic';
    data.simulation_type.numerical_method = 'fem_____';
    data.simulation_type.domain = 'freq_dom';
    data.simulation_type.dimension = '3_dimens';
    data.freq_stp = [];
    data.ni2coord = [];
    data.sei_2_ni = [];
    data.vei_2_ni = [];
    data.cpx_pres = [];
    
    % read element order
    order = strtrim(fread(fid, 8, 'char=>char')');
    if contains(order, '1st_ord_')
        element_order = 1;
        data.simulation_type.element_order = '1st_ord_';
    elseif contains(order, '2nd_ord_')
        element_order = 2;
        data.simulation_type.element_order = '2nd_ord_';
    else
        error('Unknown element order: %s', order);
    end
    
    % Read chunks until end of file
    while ~feof(fid)
        
        chunk_id = fread(fid, 8, 'char=>char')';
        if strcmp(chunk_id, 'numaveof')
            break
        end
        chunk_size = fread(fid, 1, 'uint64=>uint64');
    
        switch chunk_id
        case 'freq_stp'
            data.freq_stp = fread(fid, chunk_size/8, 'double=>double');
            
        case 'ni2coord'
            coordinate_count = chunk_size / 8;
            if mod(coordinate_count, 3) ~= 0
                error('ni2coord data length is not a multiple of 3');
            end
            coords = fread(fid, coordinate_count, 'double=>double');
            node_count = coordinate_count / 3;
            data.ni2coord = reshape(coords, 3, node_count);
            
        case 'sei_2_ni'
            ni_count = chunk_size / 8;
            % reshape based on element order
            if element_order == 1
                nodes_per_element = 3;
            else
                nodes_per_element = 6;
            end
            if mod(ni_count, nodes_per_element) ~= 0
                error('sei_2_ni data length is not a multiple of %d.', ...
                    nodes_per_element)
            end
            
            indices = fread(fid, ni_count, 'uint64=>uint64');
            element_count = ni_count / nodes_per_element;
            data.sei_2_ni = ...
                reshape(indices, nodes_per_element, element_count);
            
        case 'vei_2_ni'
            ni_count = chunk_size / 8;
            % reshape based on element order
            if element_order == 1
                nodes_per_element = 4;
            else
                nodes_per_element = 10;
            end
            if mod(ni_count, nodes_per_element) ~= 0
                error('vei_2_ni data length is not a multiple of %d.', ...
                    nodes_per_element)
            end
            
            indices = fread(fid, ni_count, 'uint64=>uint64');
            element_count = ni_count / nodes_per_element;
            data.vei_2_ni = ...
                reshape(indices, nodes_per_element, element_count);
        
        case 'cpx_pres'
            complex_count = chunk_size / 16;
            realParts = fread(fid, complex_count, 'double=>double', 8);
            fseek(fid, -int64(chunk_size)+8, 'cof');
            imagParts = fread(fid, complex_count, 'double=>double', 8);
            fseek(fid, -8, 'cof');
            data.cpx_pres = complex(realParts, imagParts);

        otherwise
            % unknown chunk type - skip it
            warning('Unknown chunk ID: %s. Skipping %d bytes.', ...
                chunk_id, chunk_size);
            fseek(fid, chunk_size, 'cof');
        end
    end
    
    % validate data consistency
    if isempty(data.freq_stp)
        error('Could not determine frequency steps.');
    end
    if isempty(data.ni2coord)
        error('Could not determine nodes.');
    end
    if isempty(data.sei_2_ni)
        error('Could not determine surface elements.');
    end
    if isempty(data.vei_2_ni)
        error('Could not determine volume elements.');
    end
    if isempty(data.cpx_pres)
        error('Could not determine pressure solution.');
    end

    % reshape to complex pressure to (ni_count × freq_count)
    node_count = size(data.ni2coord, 2);
    freq_count = length(data.freq_stp);
    if node_count * freq_count ~= length(data.cpx_pres)
        error(['cpx_pres dimensions mismatch. ' ...
            'Expected %d×%d=%d, got %d elements.'], ...
            node_count, freq_count, node_count*freq_count, complex_count);
    end
    data.cpx_pres = reshape(data.cpx_pres, node_count, freq_count);

    % convert to one-based index
    data.sei_2_ni = data.sei_2_ni + 1;
    data.vei_2_ni = data.vei_2_ni + 1;

end

end
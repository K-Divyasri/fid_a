%% MRSI spectrum functions
function tf = is_mrsi_data(p)
    % Determine if the data is MRSI spectroscopic data
    tf = false;
    
    % Check for multiple volumes (spectral points)
    nVol = size(p.nii.img, 4);
    if nVol < 2, return; end
    
    % Method 1: Check intent code for time series or other spectroscopic codes
    if any(p.nii.hdr.intent_code == [2003, 2001, 2002])
        tf = true; return;
    end
    
    % Method 2: Check if it has ppm field (common for MRSI)
    if isfield(p, 'ppm') || isfield(p.nii, 'ppm')
        tf = true; return;
    end
    
    % Method 3: Check dimensions (many spectral points, modest spatial resolution)
    dim = p.nii.hdr.dim(2:4);
    if nVol > 64 && all(dim <= 128)
        tf = true; return;
    end
    
    % Method 4: Check temporal unit for Hz or ppm
    temporal_unit = bitand(p.nii.hdr.xyzt_units, 56);
    if any(temporal_unit == [32, 40])
        tf = true; return;
    end
    
    % Method 5: Check filename
    [~, fname] = fileparts(p.nii.hdr.file_name);
    if contains(lower(fname), {'mrsi', 'spectro', 'spec', 'mrs'})
        tf = true; return;
    end
end
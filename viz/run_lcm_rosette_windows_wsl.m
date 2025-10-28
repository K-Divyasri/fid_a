function run_lcm_rosette_windows_wsl(fileName, ftSpec_smooth, ftSpec_smooth_w, brain_area_raw, figurefoldername)
% RUN_LCM_ROSETTE_WINDOWS_WSL
%  - Primary goal: produce *perfect* LCModel inputs using io_CSIwritelcm at the CSI level.
%  - Optional: also run LCModel via WSL for each brain-mask voxel using a robust bash script.
%
% Inputs:
%   fileName         : full path to the WS dataset (used only to anchor output location)
%   ftSpec_smooth    : WS CSI object (B0 corrected + masked + apodized)
%   ftSpec_smooth_w  : WU CSI object (B0 corrected + apodized)
%   brain_area_raw   : binary mask (y,x) of voxels to process
%   figurefoldername : not used here; kept for API compatibility

    if nargin < 5, figurefoldername = ''; end

    %---------------------------------------------------------------------%
    % 0) Work in the folder that holds fileName; create a clean LCM folder
    %---------------------------------------------------------------------%
    [workDir, baseName, ~] = fileparts(fileName);
    if isempty(workDir), workDir = pwd; end
    cd(workDir);

    lcmRoot = fullfile(workDir, [baseName '_lcm']);
    if ~exist(lcmRoot, 'dir'), mkdir(lcmRoot); end

    %---------------------------------------------------------------------%
    % 1) Produce LCModel inputs ONCE at the CSI level (recommended)
    %    io_CSIwritelcm will:
    %      - extract per-voxel time-domain FIDs correctly
    %      - set NUNFIL / DELTAT / HZPPPM
    %      - write valid .RAW files with correct channel ordering
    %      - emit the folder structure LCModel expects for CSI batches
    %---------------------------------------------------------------------%
    wsCSIstub = fullfile(lcmRoot, [baseName '_ftSpec_smooth_lcm']);
    wuCSIstub = fullfile(lcmRoot, [baseName '_ftSpec_smooth_lcm_w']);

    fprintf('Writing CSI-level LCModel inputs (WS/WU) under:\n  %s\n', lcmRoot);
    % These two calls are the “gold path”.
    io_CSIwritelcm(ftSpec_smooth,  wsCSIstub);
    io_CSIwritelcm(ftSpec_smooth_w, wuCSIstub);

    % NOTE: At this point, all voxel RAWs + LCModel batch layout are created.
    %       You can stop here and run LCModel separately if you want.

    %---------------------------------------------------------------------%
    % 2) OPTIONAL: also run LCModel right now via WSL for each voxel
    %---------------------------------------------------------------------%
    doRunLCM = true;   % <- set false if you only want to generate inputs

    if ~doRunLCM
        fprintf('✓ LCModel inputs ready. Skipping LCModel run (doRunLCM=false).\n');
        return;
    end

    % Your WSL script (the *fixed* one below)
    run_script_location_win = 'C:\Users\divya\Downloads\linux-send\OneDrive_1_7-5-2025\sneha3T_24x24_TE23.sh';
    run_script_location_wsl = win2wsl(run_script_location_win);

    % Brain-mask voxels (y,x)
    [loc_y, loc_x] = find(brain_area_raw == 1);
    nvox = numel(loc_x);
    fprintf('Running LCModel for %d voxels in mask...\n', nvox);

    % Result collector (single out folder)
    outCollector = fullfile(lcmRoot, [baseName '_out']);
    if ~exist(outCollector, 'dir'), mkdir(outCollector); end

    % Walk each voxel and run LCModel via WSL
    for k = 1:nvox
        vx = loc_x(k);  vy = loc_y(k);
        tag = sprintf('%dx%d_ftSpec_smooth_lcm', vx, vy);

        % Where io_CSIwritelcm put the voxel stubs (without extension)
        raw_ws = fullfile(lcmRoot, tag);
        raw_wu = [raw_ws '_w'];

        % Compute timing from the actual (per-voxel) time-domain object
        % by extracting the underlying time-domain voxel from ftSpec_smooth
        % (io_CSIwritelcm would have used these internally as well)
        sv_ws = op_CSItoMRS(ftSpec_smooth, vx, vy);
        if isfield(sv_ws,'spectralDwellTime') && ~isempty(sv_ws.spectralDwellTime)
            deltat = sv_ws.spectralDwellTime;
        else
            deltat = 1 / sv_ws.spectralWidth;
        end
        if isfield(sv_ws,'txfrq') && ~isempty(sv_ws.txfrq)
            hzpppm = sv_ws.txfrq / 1e6;    % MHz
        else
            error('Missing txfrq in voxel (%d,%d).', vx, vy);
        end
        % Time points in the FID (NUNFIL)
        if isfield(sv_ws,'dims') && isfield(sv_ws,'sz') && sv_ws.dims.t > 0
            nunfil = sv_ws.sz(sv_ws.dims.t);
        else
            nunfil = size(sv_ws.data, 1);
        end

        % Build WSL command (QUOTED paths; no unquoted $vars in the script)
        cmd = sprintf('wsl bash "%s" "%s" "%s" "%.9f" "%.6f" "%d" "%d" "%d" "2.3"', ...
            run_script_location_wsl, ...
            win2wsl(raw_ws), win2wsl(raw_wu), ...
            deltat, hzpppm, nunfil, ...
            ftSpec_smooth.sz(ftSpec_smooth.dims.x), ...
            ftSpec_smooth.sz(ftSpec_smooth.dims.y));

        fprintf('\n=== Voxel %d/%d  (%d,%d) ===\n', k, nvox, vx, vy);
        fprintf('Running LCModel:\n  %s\n', cmd);
        status = system(cmd);
        if status ~= 0
            warning('LCModel script failed for voxel (%d,%d).', vx, vy);
        end

        % Collect outputs if present
        oldlocation = [raw_ws '_out'];
        fileExt = {'.table', '.coord', '.control', '.ps', '.print', '.csv'};
        for e = 1:numel(fileExt)
            src = fullfile(oldlocation, [tag fileExt{e}]);
            if exist(src, 'file')
                fprintf('  → %s\n', src);
                movefile(src, outCollector);
            end
        end
        if exist(oldlocation, 'dir'), rmdir(oldlocation, 's'); end
    end

    fprintf('\n✔ All done. Collected LCModel outputs in:\n   %s\n', outCollector);
end

% Helper: Windows → WSL path
function linuxPath = win2wsl(winPath)
    linuxPath = strrep(winPath, '\', '/');
    linuxPath = ['/mnt/' lower(linuxPath(1)) linuxPath(3:end)];
end

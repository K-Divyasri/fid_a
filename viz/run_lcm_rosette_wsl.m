

function run_lcm_rosette_windows_wsl(fileName, ftSpec_smooth, ...
                                     ftSpec_smooth_w, brain_area_raw, ...
                                     figurefoldername)

%-------------------------------------------------------------------------%
% 0.  Where are we working?  → same folder that holds <fileName>
%-------------------------------------------------------------------------%
    [workDir, baseName, ~] = fileparts(fileName);
    if isempty(workDir), workDir = pwd; end          % if caller passed only a bare file name
    cd(workDir);

%-------------------------------------------------------------------------%
% 1.  Toolbox availability
%-------------------------------------------------------------------------%
    % addpath('C:\Users\divya\Downloads');   % holds op_CSILCModelMaps

%-------------------------------------------------------------------------%
% 2.  LCModel bash script, converted to a WSL path
%-------------------------------------------------------------------------%
    run_script_location_win = ...
        'C:\Users\divya\Downloads\linux-send\OneDrive_1_7-5-2025\sneha3T_24x24_TE2_3.sh';
    run_script_location_wsl = win2wsl(run_script_location_win);

%-------------------------------------------------------------------------%
% 3.  Find voxels we need to process
%-------------------------------------------------------------------------%
    [loc_y, loc_x] = find(brain_area_raw == 1);
    nvox = numel(loc_x);
    fprintf('Processing %d voxels in folder: %s\n', nvox, workDir);

%-------------------------------------------------------------------------%
% 4.  Loop over voxels
%-------------------------------------------------------------------------%
    for k = 1:nvox
        vx = loc_x(k);   vy = loc_y(k);
        tag = sprintf('%dx%d_ftSpec_smooth_lcm', vx, vy);

        fprintf('\n=== Voxel %d of %d  (%d,%d) ===\n', k, nvox, vx, vy);

    %--- 4a.  Write .RAW files directly in workDir -----------------------%
        rawfile1 = fullfile(workDir, tag);
        rawfile2 = [rawfile1 '_w'];

        fprintf('Writing .RAW files:\n  %s\n  %s\n', rawfile1, rawfile2);
        io_writelcm(op_CSItoMRS(ftSpec_smooth,   vx, vy), rawfile1, ftSpec_smooth.te);
        io_writelcm(op_CSItoMRS(ftSpec_smooth_w, vx, vy), rawfile2, ftSpec_smooth_w.te);

    %--- 4b.  Call LCModel via WSL --------------------------------------%
        cmd = sprintf('wsl bash %s %s %s', ...
            run_script_location_wsl, win2wsl(rawfile1), win2wsl(rawfile2));
        fprintf('Running LCModel:\n  %s\n', cmd);
        if system(cmd) ~= 0
            warning('LCModel script failed for voxel (%d,%d).', vx, vy);
        end

    %--- 4c.  Move LCModel result files back to workDir ------------------%
        oldlocation = [rawfile1 '_out'];         % LCModel’s auto folder
        fileExt = {'.table', '.coord', '.control', '.ps', '.print'};

        for e = 1:numel(fileExt)
            src = fullfile(oldlocation, [tag fileExt{e}]);
            if exist(src, 'file')
                fprintf('Moving %s → %s\n', src, workDir);
                movefile(src, workDir);
            else
                warning('Expected file not found: %s', src);
            end
        end

        if exist(oldlocation, 'dir')
            rmdir(oldlocation, 's');
        end
    end
end
%-------------------------------------------------------------------------%
% 5.  Build metabolite maps
%-------------------------------------------------------------------------%
% fprintf('\nAll voxels processed.  Generating maps from: %s\n', workDir);
% addpath(genpath('C:\Users\divya\Downloads\linux-send\fida2\new'));
% 
% final_location = 'C:\Users\divya\Downloads\Lexar\DataForDivya\meas_MID01415_FID95676_lb_RosetteSpinEcho_ek_v1a_lcm\meas_MID01415_FID95676_lb_RosetteSpinEcho_ek_v1a_out';
% figurefoldername = 'Figures_meas_MID01415_FID95676_lb_RosetteSpinEcho_ek_v1a';
% 
% [map, crlb, LW, SNR] = op_CSILCModelMaps(40,40,final_location,'figure_folder_name',figurefoldername);
% 
% % Visualise (example)
% names = {'CrPCr','GPCPCh','NAANAAG','GluGln','Ins'};
% for i=1:numel(names)
%     figure('Name',names{i});
%     imagesc(map.(names{i})); axis image; colorbar; colormap hot;
%     title(['LCModel ' names{i}]);
% end
%     if exist('op_CSILCModelMaps', 'file')
%         [map, crlb, LW, SNR] = op_CSILCModelMaps(48, 48, workDir, ...
%                                    'figure_folder_name', figurefoldername);
%     else
%         warning('op_CSILCModelMaps is not on the MATLAB path.');
%     end
% end

%--- Helper: Windows -> WSL path -----------------------------------------%
function linuxPath = win2wsl(winPath)
    linuxPath = strrep(winPath, '\', '/');
    linuxPath = ['/mnt/' lower(linuxPath(1)) linuxPath(3:end)];
end

% function run_lcm_rosette_windows_wsl(fileName, ftSpec_smooth, ftSpec_smooth_w, brain_area_raw, figurefoldername)
%     % --- Setup paths ---
%     addpath('C:\Users\divya\Downloads'); % Ensure op_CSILCModelMaps is available
%     rootdir = pwd;
% 
%     % --- LCModel WSL script ---
%     run_script_location_win = 'C:\Users\divya\Downloads\fida2\new\sneha3T_48x48_TE20.sh';
%     run_script_location_wsl = win2wsl(run_script_location_win);
% 
%     % --- Prepare LCM output directories ---
%     [~, baseName, ~] = fileparts(fileName);
%     lcmfoldername = fullfile(rootdir, [baseName '_lcm']);
%     if ~exist(lcmfoldername, 'dir'), mkdir(lcmfoldername); end
% 
%     lcmoutput = [baseName '_out'];
%     outfolder = fullfile(lcmfoldername, lcmoutput);
%     if ~exist(outfolder, 'dir'), mkdir(outfolder); end
% 
%     % --- Locate voxels ---
%     [loc_y, loc_x] = find(brain_area_raw==1);
%     fprintf('Processing %d voxels...\n', numel(loc_x));
% 
%     for k = 1:numel(loc_x)
%         fprintf('\n=== Voxel %d (%d, %d) ===\n', k, loc_x(k), loc_y(k));
%         lcm_svs_file_name = sprintf('%dx%d_ftSpec_smooth_lcm', loc_x(k), loc_y(k));
%         lcm_svs_file_name_w = [lcm_svs_file_name '_w'];
% 
%         % --- Write .RAW files ---
%         rawfile1 = fullfile(lcmfoldername, lcm_svs_file_name);
%         rawfile2 = fullfile(lcmfoldername, lcm_svs_file_name_w);
%         fprintf('Writing .RAW files: %s, %s\n', rawfile1, rawfile2);
%         io_writelcm(op_CSItoMRS(ftSpec_smooth,   loc_x(k), loc_y(k)), rawfile1, ftSpec_smooth.te);
%         io_writelcm(op_CSItoMRS(ftSpec_smooth_w, loc_x(k), loc_y(k)), rawfile2, ftSpec_smooth_w.te);
% 
%         % --- Run LCModel ---
%         vox_lcm   = win2wsl(rawfile1);
%         vox_lcm_w = win2wsl(rawfile2);
%         command = sprintf('wsl bash %s %s %s', run_script_location_wsl, vox_lcm, vox_lcm_w);
%         fprintf('Running LCModel in WSL:\n  %s\n', command);
%         status = system(command);
%         if status ~= 0
%             warning('LCModel script failed for voxel %d (%d,%d)', k, loc_x(k), loc_y(k));
%         end
% 
%         % --- Move LCModel output files ---
%         oldlocation = fullfile(lcmfoldername, [lcm_svs_file_name '_out']);
%         variable_type = {'.table', '.coord', '.control', '.ps'};
% 
%         for i = 1:numel(variable_type)
%             file_to_move = fullfile(oldlocation, [lcm_svs_file_name variable_type{i}]);
%             if exist(file_to_move, 'file')
%                 fprintf('Moving file: %s → %s\n', file_to_move, outfolder);
%                 movefile(file_to_move, outfolder);
%             else
%                 warning('Could not find file: %s', file_to_move);
%             end
%         end
% 
%         if exist(oldlocation, 'dir')
%             fprintf('Removing temporary folder: %s\n', oldlocation);
%             rmdir(oldlocation, 's');
%         end
%     end
% 
%     % --- Generate maps ---
%     final_location = outfolder;
%     fprintf('\nAll voxels processed. Generating LCModel maps from: %s\n', final_location);
% 
%     if exist('op_CSILCModelMaps', 'file')
%         [map, crlb, LW, SNR] = op_CSILCModelMaps(48, 48, final_location, 'figure_folder_name', figurefoldername);
%     else
%         warning('op_CSILCModelMaps is not on the MATLAB path.');
%     end
% end
% 
% function linuxPath = win2wsl(winPath)
%     linuxPath = strrep(winPath, '\', '/');
%     drive = lower(linuxPath(1));
%     linuxPath = ['/mnt/' drive linuxPath(3:end)];
% end


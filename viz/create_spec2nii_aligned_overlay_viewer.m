function create_spec2nii_aligned_overlay_viewer(ftSpec, map, crlb, LW, SNR, t1_path, spec2nii_empty_path)
% MRSI viewer with CORRECTED orientation matching spec2nii
    
    fprintf('\n========================================================================\n');
    fprintf('=== CREATING SPEC2NII ALIGNED MRSI OVERLAY ===\n');
    fprintf('========================================================================\n\n');
    
    % Load T1
    t1 = nii_tool('load', t1_path);
    metabolite_names = fieldnames(map);
    
    % Get T1 spatial transformation
    if t1.hdr.sform_code > 0
        if isfield(t1.hdr, 'sform_mat')
            T1_affine = [double(t1.hdr.sform_mat); 0 0 0 1];
        else
            T1_affine = [double(t1.hdr.srow_x); 
                        double(t1.hdr.srow_y); 
                        double(t1.hdr.srow_z); 
                        0 0 0 1];
        end
    else
        error('T1 must have spatial transform (sform)');
    end
    
    % Load spec2nii empty NIfTI
    fprintf('Loading spec2nii empty NIfTI...\n');
    if ~exist(spec2nii_empty_path, 'file')
        error('spec2nii empty NIfTI not found: %s', spec2nii_empty_path);
    end
    
    nii_spec2nii = nii_tool('hdr', spec2nii_empty_path);
    
    % Extract spec2nii positioning
    if nii_spec2nii.sform_code > 0
        if isfield(nii_spec2nii, 'sform_mat') && ~isempty(nii_spec2nii.sform_mat)
            spec2nii_sform = [double(nii_spec2nii.sform_mat); 0 0 0 1];
        else
            error('Could not extract sform from spec2nii NIfTI');
        end
        
        % Calculate spec2nii center
        spec2nii_dims = double(nii_spec2nii.dim(2:4));
        spec2nii_center_vox = (spec2nii_dims - 1) / 2;
        spec2nii_center_world = spec2nii_sform * [spec2nii_center_vox'; 1];
        spec2nii_slice_center = spec2nii_center_world(1:3);
        
        % Extract orientation vectors
        x_vec_spec2nii = spec2nii_sform(1:3,1);
        y_vec_spec2nii = spec2nii_sform(1:3,2);
        z_vec_spec2nii = spec2nii_sform(1:3,3);
        
        x_dir_spec2nii = x_vec_spec2nii / norm(x_vec_spec2nii);
        y_dir_spec2nii = y_vec_spec2nii / norm(y_vec_spec2nii);
        z_dir_spec2nii = z_vec_spec2nii / norm(z_vec_spec2nii);
        
    else
        error('spec2nii NIfTI has no spatial transform');
    end
    
    % Get ftSpec information
    [mr, mc] = size(map.(metabolite_names{1}));
    
    % Extract ftSpec vectors
    ftspec_x_vec = ftSpec.affineMatrix(1:3,1);
    ftspec_y_vec = ftSpec.affineMatrix(1:3,2);
    ftspec_z_vec = ftSpec.affineMatrix(1:3,3);
    
    dx_ftspec = norm(ftspec_x_vec);
    dy_ftspec = norm(ftspec_y_vec);
    dz_ftspec = norm(ftspec_z_vec);
    
    ftspec_x_dir = ftspec_x_vec / dx_ftspec;
    ftspec_y_dir = ftspec_y_vec / dy_ftspec;
    
    fprintf('MRSI grid: %d x %d voxels\n', mr, mc);
    fprintf('MRSI voxel sizes: [%.2f, %.2f, %.2f] mm\n', dx_ftspec, dy_ftspec, dz_ftspec);
    
    % Check orientation alignment
    dot_x = dot(ftspec_x_dir, x_dir_spec2nii);
    dot_y = dot(ftspec_y_dir, y_dir_spec2nii);
    
    fprintf('\n=== ORIENTATION MATCHING ===\n');
    fprintf('ftSpec X: [%.3f, %.3f, %.3f]\n', ftspec_x_dir(1), ftspec_x_dir(2), ftspec_x_dir(3));
    fprintf('spec2nii X: [%.3f, %.3f, %.3f]\n', x_dir_spec2nii(1), x_dir_spec2nii(2), x_dir_spec2nii(3));
    fprintf('Dot product X: %.3f %s\n', dot_x, ternary(dot_x < 0, '(OPPOSITE - will flip)', '(same)'));
    fprintf('ftSpec Y: [%.3f, %.3f, %.3f]\n', ftspec_y_dir(1), ftspec_y_dir(2), ftspec_y_dir(3));
    fprintf('spec2nii Y: [%.3f, %.3f, %.3f]\n', y_dir_spec2nii(1), y_dir_spec2nii(2), y_dir_spec2nii(3));
    fprintf('Dot product Y: %.3f %s\n', dot_y, ternary(dot_y < 0, '(OPPOSITE - will flip)', '(same)'));
    
    % CRITICAL FIX: Use spec2nii orientation with ftSpec voxel sizes
    % This ensures the overlay has the same orientation as spec2nii
    target_center = spec2nii_slice_center;
    
    % Build vectors using spec2nii directions and ftSpec magnitudes
    x_vec_final = x_dir_spec2nii * dx_ftspec;
    y_vec_final = y_dir_spec2nii * dy_ftspec;
    z_vec_final = z_dir_spec2nii * dz_ftspec;
    
    fprintf('\n=== FINAL VECTORS (spec2nii orientation + ftSpec sizes) ===\n');
    fprintf('X vector: [%.4f, %.4f, %.4f], mag=%.2f\n', x_vec_final(1), x_vec_final(2), x_vec_final(3), norm(x_vec_final));
    fprintf('Y vector: [%.4f, %.4f, %.4f], mag=%.2f\n', y_vec_final(1), y_vec_final(2), y_vec_final(3), norm(y_vec_final));
    fprintf('Z vector: [%.4f, %.4f, %.4f], mag=%.2f\n', z_vec_final(1), z_vec_final(2), z_vec_final(3), norm(z_vec_final));
    
    % Calculate corner position
    half_grid_x = (mc-1)/2 * x_vec_final;
    half_grid_y = (mr-1)/2 * y_vec_final;
    corner_pos = target_center - half_grid_x - half_grid_y;
    
    fprintf('\n=== POSITIONING ===\n');
    fprintf('Target center: [%.2f, %.2f, %.2f]\n', target_center(1), target_center(2), target_center(3));
    fprintf('Corner position: [%.2f, %.2f, %.2f]\n', corner_pos(1), corner_pos(2), corner_pos(3));
    
    % Build final affine
    MRSI_affine = eye(4);
    MRSI_affine(1:3,1) = x_vec_final;
    MRSI_affine(1:3,2) = y_vec_final;
    MRSI_affine(1:3,3) = z_vec_final;
    MRSI_affine(1:3,4) = corner_pos;
    
    % Verify
    verify_center = corner_pos + half_grid_x + half_grid_y;
    fprintf('Verification error: %.6f mm\n', norm(verify_center - target_center));
    
    fprintf('\n=== TEST CORNER POSITIONS ===\n');
    test_corners = [0, 0; mc-1, 0; mc-1, mr-1; 0, mr-1];
    for k = 1:4
        test_pos = [test_corners(k,:), 0, 1]';
        world_pos = MRSI_affine * test_pos;
        fprintf('Voxel [%2d,%2d] -> world: [%7.2f, %7.2f, %6.2f]\n', ...
            test_corners(k,1), test_corners(k,2), world_pos(1), world_pos(2), world_pos(3));
    end
    
    % Find optimal slice
    center_t1_vox = inv(T1_affine) * [target_center; 1];
    optimal_slice = round(center_t1_vox(3) + 1);
    fprintf('\nOptimal T1 slice: %d\n', optimal_slice);
    fprintf('========================================================================\n\n');
    
    % Create figure
    fig = figure('Name', 'spec2nii Aligned MRSI Overlay', 'Position', [100 100 1400 900], 'Color', 'k');
    
    current_slice = optimal_slice;
    current_met = 1;
    
    % Controls
    uicontrol('Style', 'text', 'String', 'Metabolite:', ...
              'Position', [20 850 80 25], 'BackgroundColor', 'k', 'ForegroundColor', 'w');
    met_popup = uicontrol('Style', 'popupmenu', 'String', metabolite_names, ...
                          'Position', [110 850 120 25], 'Callback', @update_view);
    
    uicontrol('Style', 'text', 'String', 'Slice:', ...
              'Position', [250 850 50 25], 'BackgroundColor', 'k', 'ForegroundColor', 'w');
    slice_slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', size(t1.img, 3), ...
                            'Value', current_slice, 'Position', [310 850 200 25], ...
                            'Callback', @update_view);
    slice_text = uicontrol('Style', 'text', 'String', sprintf('%d/%d', current_slice, size(t1.img,3)), ...
                          'Position', [520 850 50 25], 'BackgroundColor', 'k', 'ForegroundColor', 'w');
    
    uicontrol('Style', 'text', 'String', 'Scaling:', ...
              'Position', [600 850 60 25], 'BackgroundColor', 'k', 'ForegroundColor', 'w');
    scale_popup = uicontrol('Style', 'popupmenu', ...
                           'String', {'Auto (95th %ile)', 'Max value', 'Median x2'}, ...
                           'Position', [670 850 130 25], 'Callback', @update_view);
    
    uicontrol('Style', 'text', 'String', 'Slice tol:', ...
              'Position', [820 850 70 25], 'BackgroundColor', 'k', 'ForegroundColor', 'w');
    tol_edit = uicontrol('Style', 'edit', 'String', '10', ...
                         'Position', [900 850 40 25], 'Callback', @update_view);
    
    ax = axes('Position', [0.1 0.1 0.65 0.8]);
    
    info_text = uicontrol('Style', 'text', ...
                         'Position', [950 100 400 750], ...
                         'BackgroundColor', [0.1 0.1 0.1], ...
                         'ForegroundColor', 'w', ...
                         'FontName', 'FixedWidth', ...
                         'FontSize', 10, ...
                         'HorizontalAlignment', 'left', ...
                         'String', 'Click on MRSI voxel');
    
    % Store data
    data.t1 = t1;
    data.map = map;
    data.crlb = crlb;
    data.LW = LW;
    data.SNR = SNR;
    data.metabolite_names = metabolite_names;
    data.ax = ax;
    data.met_popup = met_popup;
    data.slice_slider = slice_slider;
    data.slice_text = slice_text;
    data.scale_popup = scale_popup;
    data.tol_edit = tol_edit;
    data.info_text = info_text;
    data.T1_affine = T1_affine;
    data.T1_inv = inv(T1_affine);
    data.MRSI_affine = MRSI_affine;
    data.mr = mr;
    data.mc = mc;
    
    set(fig, 'UserData', data);
    set(fig, 'WindowButtonDownFcn', @click_handler);
    
    update_view();
    
    function update_view(~, ~)
        data = get(fig, 'UserData');
        current_met = get(data.met_popup, 'Value');
        current_slice = round(get(data.slice_slider, 'Value'));
        scale_mode = get(data.scale_popup, 'Value');
        slice_tolerance = str2double(get(data.tol_edit, 'String'));
        
        set(data.slice_text, 'String', sprintf('%d/%d', current_slice, size(data.t1.img,3)));
        
        % Get T1 slice
        t1_slice = double(squeeze(data.t1.img(:,:,current_slice)));
        t1_slice = (t1_slice - min(t1_slice(:))) / (max(t1_slice(:)) - min(t1_slice(:)));
        
        % Get metabolite data
        met_name = data.metabolite_names{current_met};
        met_data = data.map.(met_name);
        
        % Determine color scaling
        nonzero_vals = met_data(met_data > 0);
        if isempty(nonzero_vals)
            max_val = 1e-12;
        else
            switch scale_mode
                case 1
                    max_val = prctile(nonzero_vals, 95);
                case 2
                    max_val = max(nonzero_vals);
                case 3
                    max_val = median(nonzero_vals) * 2;
            end
        end
        
        if max_val == 0
            max_val = 1e-12;
        end
        
        % Display T1
        axes(data.ax);
        cla;
        imagesc(t1_slice');
        colormap(data.ax, gray);
        axis image; axis off; hold on;
        
        % Overlay MRSI voxels
        voxels_drawn = 0;
        
        for i = 1:data.mr
            for j = 1:data.mc
                val = met_data(i,j);
                
                if val > 0
                    % Voxel corners in MRSI coordinate system (0-based)
                    corners_mrsi = [
                        j-1, i-1, 0, 1;
                        j,   i-1, 0, 1;
                        j,   i,   0, 1;
                        j-1, i,   0, 1
                    ]';
                    
                    % Transform to world coordinates
                    corners_world = data.MRSI_affine * corners_mrsi;
                    
                    % Transform to T1 voxel coordinates
                    corners_t1 = data.T1_inv * corners_world;
                    
                    % Extract coordinates (convert to 1-based)
                    x_coords = corners_t1(1,:) + 1;
                    y_coords = corners_t1(2,:) + 1;
                    z_coords = corners_t1(3,:) + 1;
                    
                    % Check if voxel is near current slice
                    z_dist = abs(mean(z_coords) - current_slice);
                    
                    if z_dist < slice_tolerance
                        % Normalize intensity
                        intensity = val / max_val;
                        intensity = max(0, min(1, intensity));
                        
                        % Hot colormap
                        if intensity < 0.33
                            color = [intensity*3, 0, 0];
                        elseif intensity < 0.67
                            color = [1, (intensity-0.33)*3, 0];
                        else
                            color = [1, 1, (intensity-0.67)*3];
                        end
                        color = max(0, min(1, color));
                        
                        % Draw filled polygon
                        if all(isfinite(x_coords)) && all(isfinite(y_coords))
                            try
                                patch(x_coords, y_coords, color, ...
                                      'EdgeColor', 'w', 'LineWidth', 0.5, ...
                                      'FaceAlpha', 0.7);
                                voxels_drawn = voxels_drawn + 1;
                            catch
                                % Skip problematic voxels
                            end
                        end
                    end
                end
            end
        end
        
        title(sprintf('%s (Slice %d, tol=%.0f) - Max: %.2e - Voxels: %d', ...
                     met_name, current_slice, slice_tolerance, max_val, voxels_drawn), ...
              'Color', 'w', 'FontSize', 12);
        
        set(fig, 'UserData', data);
    end
    
    function click_handler(~, ~)
        data = get(fig, 'UserData');
        
        pos = get(data.ax, 'CurrentPoint');
        x_t1 = pos(1,1);
        y_t1 = pos(1,2);
        z_t1 = round(get(data.slice_slider, 'Value'));
        
        % Convert to world then MRSI coordinates
        world_pos = data.T1_affine * [x_t1-1; y_t1-1; z_t1-1; 1];
        mrsi_pos = inv(data.MRSI_affine) * world_pos;
        
        col = round(mrsi_pos(1)) + 1;
        row = round(mrsi_pos(2)) + 1;
        
        fprintf('\nClick: T1[%.1f,%.1f,%d] -> World[%.2f,%.2f,%.2f] -> MRSI(%d,%d)\n', ...
            x_t1, y_t1, z_t1, world_pos(1), world_pos(2), world_pos(3), row, col);
        
        if row >= 1 && row <= data.mr && col >= 1 && col <= data.mc
            info_str = sprintf('VOXEL (%d, %d)\n%s\n\n', row, col, repmat('=', 1, 30));
            info_str = [info_str 'METABOLITE CONCENTRATIONS:\n'];
            
            for m = 1:length(data.metabolite_names)
                met_name = data.metabolite_names{m};
                conc = data.map.(met_name)(row, col);
                crlb_val = data.crlb.(met_name)(row, col);
                
                if crlb_val <= 20
                    quality = 'Good';
                elseif crlb_val <= 50
                    quality = 'Fair';
                else
                    quality = 'Poor';
                end
                
                info_str = [info_str sprintf('%s: %.3e\n  CRLB: %.1f%% (%s)\n', ...
                           met_name, conc, crlb_val, quality)];
            end
            
            lw_val = data.LW(row, col);
            snr_val = data.SNR(row, col);
            
            if lw_val <= 0.1
                lw_quality = 'Good';
            else
                lw_quality = 'Poor';
            end
            
            if snr_val >= 5
                snr_quality = 'Good';
            else
                snr_quality = 'Poor';
            end
            
            info_str = [info_str newline 'QUALITY METRICS:' newline];
            info_str = [info_str sprintf('Linewidth: %.4f ppm (%s)\n', lw_val, lw_quality)];
            info_str = [info_str sprintf('SNR: %.2f (%s)', snr_val, snr_quality)];
            
            set(data.info_text, 'String', info_str);
            
            % Highlight clicked voxel
            axes(data.ax);
            delete(findobj(data.ax, 'Tag', 'highlight'));
            
            corners_mrsi = [col-1, row-1, 0, 1; col, row-1, 0, 1; col, row, 0, 1; col-1, row, 0, 1]';
            corners_world = data.MRSI_affine * corners_mrsi;
            corners_t1 = data.T1_inv * corners_world;
            
            x_coords = corners_t1(1,:) + 1;
            y_coords = corners_t1(2,:) + 1;
            
            if all(isfinite(x_coords)) && all(isfinite(y_coords))
                plot([x_coords, x_coords(1)], [y_coords, y_coords(1)], ...
                     'c-', 'LineWidth', 3, 'Tag', 'highlight');
            end
        else
            fprintf('Click outside MRSI grid\n');
        end
    end
end

function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end


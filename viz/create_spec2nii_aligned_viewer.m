function create_spec2nii_aligned_viewer(ftSpec, map, crlb, LW, SNR, t1_path, spec2nii_empty_path, y_offset_mm)
    % MRSI viewer with adjustable Y offset to shift metabolite map
        
        if nargin < 8
            y_offset_mm = 10;  % Default: shift down by 40mm
        end
        
        fprintf('\n========================================================================\n');
        fprintf('=== CREATING SPEC2NII ALIGNED MRSI OVERLAY ===\n');
        fprintf('=== T1 ROTATED 180° | Y OFFSET: %.1f mm ===\n', y_offset_mm);
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
        
        fprintf('T1 dimensions: [%d x %d x %d]\n', size(t1.img,1), size(t1.img,2), size(t1.img,3));
        
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
        
        fprintf('MRSI grid: %d x %d voxels\n', mr, mc);
        fprintf('MRSI voxel sizes: [%.2f, %.2f, %.2f] mm\n', dx_ftspec, dy_ftspec, dz_ftspec);
        
        % APPLY Y OFFSET to shift metabolite map
        target_center = spec2nii_slice_center;
        target_center(2) = target_center(2) + y_offset_mm;  % Shift in Y direction
        
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
        
        fprintf('\n=== POSITIONING (WITH Y OFFSET) ===\n');
        fprintf('Original center: [%.2f, %.2f, %.2f]\n', spec2nii_slice_center(1), spec2nii_slice_center(2), spec2nii_slice_center(3));
        fprintf('Adjusted center: [%.2f, %.2f, %.2f] (Y offset: %.1f mm)\n', target_center(1), target_center(2), target_center(3), y_offset_mm);
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
        
        % Find optimal slice
        center_t1_vox = inv(T1_affine) * [target_center; 1];
        optimal_slice = round(center_t1_vox(3) + 1);
        fprintf('\nOptimal T1 slice: %d\n', optimal_slice);
        fprintf('========================================================================\n\n');
        
        % Create figure
        fig = figure('Name', sprintf('MRSI Overlay (Y offset: %.0f mm)', y_offset_mm), ...
                     'Position', [100 100 1400 900], 'Color', 'k');
        
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
        tol_slider = uicontrol('Style', 'slider', 'Min', 0.5, 'Max', 10, ...
                              'Value', 5, 'Position', [900 850 150 25], ...
                              'Callback', @update_view);
        tol_text = uicontrol('Style', 'text', 'String', '5.0', ...
                            'Position', [1060 850 40 25], 'BackgroundColor', 'k', 'ForegroundColor', 'w');
        
        % Status
        uicontrol('Style', 'text', 'String', sprintf('Y OFFSET: %.0fmm', y_offset_mm), ...
                  'Position', [1120 850 180 25], 'BackgroundColor', 'k', ...
                  'ForegroundColor', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
        
        % Main axes
        ax = axes('Position', [0.05 0.05 0.6 0.85], 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
        
        % Info panel
        info_text = uicontrol('Style', 'text', 'String', 'Click on a voxel...', ...
                             'Position', [950 50 400 750], ...
                             'BackgroundColor', 'k', 'ForegroundColor', 'w', ...
                             'HorizontalAlignment', 'left', 'FontName', 'Courier', ...
                             'FontSize', 10);
        
        % Store data
        data = struct();
        data.t1 = t1;
        data.map = map;
        data.crlb = crlb;
        data.LW = LW;
        data.SNR = SNR;
        data.metabolite_names = metabolite_names;
        data.T1_affine = T1_affine;
        data.T1_inv = inv(T1_affine);
        data.MRSI_affine = MRSI_affine;
        data.mr = mr;
        data.mc = mc;
        data.ax = ax;
        data.slice_slider = slice_slider;
        data.slice_text = slice_text;
        data.met_popup = met_popup;
        data.scale_popup = scale_popup;
        data.tol_slider = tol_slider;
        data.tol_text = tol_text;
        data.info_text = info_text;
        
        set(fig, 'UserData', data);
        set(fig, 'WindowButtonDownFcn', @click_handler);
        
        update_view();
        
        fprintf('========================================================================\n');
        fprintf('VIEWER READY - Y offset applied: %.1f mm\n', y_offset_mm);
        fprintf('========================================================================\n');
        
        function update_view(~, ~)
            data = get(fig, 'UserData');
            
            current_slice = round(get(data.slice_slider, 'Value'));
            current_met = get(data.met_popup, 'Value');
            scale_mode = get(data.scale_popup, 'Value');
            slice_tolerance = get(data.tol_slider, 'Value');
            
            set(data.slice_text, 'String', sprintf('%d/%d', current_slice, size(data.t1.img, 3)));
            set(data.tol_text, 'String', sprintf('%.1f', slice_tolerance));
            
            met_name = data.metabolite_names{current_met};
            met_data = data.map.(met_name);
            
            % Get T1 slice and rotate 180°
            t1_slice = data.t1.img(:,:,current_slice);
            t1_slice_rotated = rot90(t1_slice, 2);
            
            % Store dimensions
            t1_width = size(data.t1.img, 1);
            t1_height = size(data.t1.img, 2);
            
            % Calculate scaling
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
            
            % Display T1 (ROTATED 180°)
            axes(data.ax);
            cla;
            imagesc(t1_slice_rotated');
            colormap(data.ax, gray);
            axis image; axis off; hold on;
            
            % Overlay MRSI voxels
            voxels_drawn = 0;
            
            for i = 1:data.mr
                for j = 1:data.mc
                    val = met_data(i,j);
                    
                    if val > 0
                        corners_mrsi = [
                            j-1, i-1, 0, 1;
                            j,   i-1, 0, 1;
                            j,   i,   0, 1;
                            j-1, i,   0, 1
                        ]';
                        
                        corners_world = data.MRSI_affine * corners_mrsi;
                        corners_t1 = data.T1_inv * corners_world;
                        
                        x_coords_orig = corners_t1(1,:) + 1;
                        y_coords_orig = corners_t1(2,:) + 1;
                        z_coords = corners_t1(3,:) + 1;
                        
                        % Apply 180° rotation
                        x_coords = t1_width + 1 - x_coords_orig;
                        y_coords = t1_height + 1 - y_coords_orig;
                        
                        z_dist = abs(mean(z_coords) - current_slice);
                        
                        if z_dist < slice_tolerance
                            intensity = val / max_val;
                            intensity = max(0, min(1, intensity));
                            
                            if intensity < 0.33
                                color = [intensity*3, 0, 0];
                            elseif intensity < 0.67
                                color = [1, (intensity-0.33)*3, 0];
                            else
                                color = [1, 1, (intensity-0.67)*3];
                            end
                            color = max(0, min(1, color));
                            
                            if all(isfinite(x_coords)) && all(isfinite(y_coords))
                                try
                                    patch(x_coords, y_coords, color, ...
                                          'EdgeColor', 'w', 'LineWidth', 0.5, ...
                                          'FaceAlpha', 0.7);
                                    voxels_drawn = voxels_drawn + 1;
                                catch
                                end
                            end
                        end
                    end
                end
            end
            
            title(sprintf('%s | Slice %d | Voxels: %d', ...
                         met_name, current_slice, voxels_drawn), ...
                  'Color', 'w', 'FontSize', 12);
            
            set(fig, 'UserData', data);
        end
        
        function click_handler(~, ~)
            data = get(fig, 'UserData');
            
            pos = get(data.ax, 'CurrentPoint');
            x_t1_display = pos(1,1);
            y_t1_display = pos(1,2);
            z_t1 = round(get(data.slice_slider, 'Value'));
            
            % Reverse 180° rotation
            t1_width = size(data.t1.img, 1);
            t1_height = size(data.t1.img, 2);
            x_t1_orig = t1_width + 1 - x_t1_display;
            y_t1_orig = t1_height + 1 - y_t1_display;
            
            world_pos = data.T1_affine * [x_t1_orig-1; y_t1_orig-1; z_t1-1; 1];
            mrsi_pos = inv(data.MRSI_affine) * world_pos;
            
            col = round(mrsi_pos(1)) + 1;
            row = round(mrsi_pos(2)) + 1;
            
            fprintf('\nClick: T1[%.1f,%.1f,%d] -> World[%.2f,%.2f,%.2f] -> MRSI(%d,%d)\n', ...
                x_t1_display, y_t1_display, z_t1, world_pos(1), world_pos(2), world_pos(3), row, col);
            
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
                
                axes(data.ax);
                delete(findobj(data.ax, 'Tag', 'highlight'));
                
                corners_mrsi = [col-1, row-1, 0, 1; col, row-1, 0, 1; col, row, 0, 1; col-1, row, 0, 1]';
                corners_world = data.MRSI_affine * corners_mrsi;
                corners_t1 = data.T1_inv * corners_world;
                
                x_coords_orig = corners_t1(1,:) + 1;
                y_coords_orig = corners_t1(2,:) + 1;
                
                x_coords = t1_width + 1 - x_coords_orig;
                y_coords = t1_height + 1 - y_coords_orig;
                
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
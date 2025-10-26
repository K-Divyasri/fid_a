function show_mrsi_spectrum_onclick(hs, p, c, ax)
    % Show MRSI spectrum in the 4th quadrant of nii_viewer
    % LEFT: Time domain (ccav_w) - real part
    % RIGHT: Frequency domain - real part
    
    fprintf('=== Starting spectrum display ===\n');
    
    % Get clicked voxel coordinates (these are in background image space)
    ijk_bg = cell2mat(get(hs.ijk, 'Value'));
    fprintf('Background voxel coordinates: [%d, %d, %d]\n', ijk_bg(1), ijk_bg(2), ijk_bg(3));
    
    try
        % Transform coordinates to MRSI image space
        if isfield(p, 'R0') % overlay with different resolution
            xyz = hs.bg.R * [ijk_bg-1; 1];
            xyz = xyz(1:3);
            fprintf('XYZ coordinates: [%.2f, %.2f, %.2f]\n', xyz(1), xyz(2), xyz(3));
            
            ijk_mrsi = round(p.Ri * [xyz; 1]) + 1;
            ijk_mrsi = ijk_mrsi(1:3);
        else
            ijk_mrsi = ijk_bg;
        end
        
        fprintf('MRSI voxel coordinates: [%d, %d, %d]\n', ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3));
        
        % Check bounds
        img_size = size(p.nii.img);
        fprintf('MRSI image size: %dx%dx%dx%d\n', img_size(1), img_size(2), img_size(3), img_size(4));
        
        if any(ijk_mrsi < 1) || ijk_mrsi(1) > img_size(1) || ijk_mrsi(2) > img_size(2) || ijk_mrsi(3) > img_size(3)
            fprintf('Coordinates out of bounds, skipping\n');
            return;
        end
        
        % Check if we have the ftSpec_smooth_w data stored
        if ~isfield(p, 'ftSpec_smooth_w') || isempty(p.ftSpec_smooth_w)
            fprintf('WARNING: ftSpec_smooth_w not found. Will only show time domain.\n');
            show_freq = false;
        else
            show_freq = true;
        end
        
        fprintf('Extracting time domain data (ccav_w)...\n');
        fid = squeeze(p.nii.img(ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3), :));
        fid = fid(:); % Ensure column vector
        
        % Check if we got valid data
        if all(fid == 0) || all(isnan(fid))
            fprintf('No signal at this voxel, skipping\n');
            return;
        end
        
        % Apply scaling if needed
        if isfield(p, 'scl_slope') && p.scl_slope ~= 0
            fid = single(fid) * p.scl_slope + p.scl_inter;
        end
        
        % Get time axis
        if isfield(p.nii.hdr, 'pixdim') && p.nii.hdr.pixdim(5) > 0
            dwelltime = double(p.nii.hdr.pixdim(5));
        else
            dwelltime = 5e-6;
        end
        
        sz = length(fid);
        t = (0:sz-1)' * dwelltime;
        t_ms = t * 1000; % Convert to milliseconds
        
        % Extract frequency domain voxel using op_CSItoMRS if available
        if show_freq
            fprintf('Extracting frequency domain voxel using op_CSItoMRS...\n');
            try
                vox_MRS = op_CSItoMRS(p.ftSpec_smooth_w, ijk_mrsi(1), ijk_mrsi(2));
                fprintf('Successfully extracted voxel using op_CSItoMRS\n');
            catch ME
                fprintf('ERROR in op_CSItoMRS: %s\n', ME.message);
                show_freq = false;
            end
        end
        
        fprintf('Creating spectrum display in 4th quadrant...\n');
        
        if ishandle(hs.ax(4))
            % Store original position if not already stored
            if ~isappdata(hs.ax(4), 'OriginalPosition')
                origPos = get(hs.ax(4), 'Position');
                setappdata(hs.ax(4), 'OriginalPosition', origPos);
            end
            
            % Restore original position
            origPos = getappdata(hs.ax(4), 'OriginalPosition');
            
            % Delete any existing subplots
            if isappdata(hs.ax(4), 'TimeAxes')
                old_ax = getappdata(hs.ax(4), 'TimeAxes');
                if ishandle(old_ax), delete(old_ax); end
            end
            if isappdata(hs.ax(4), 'FreqAxes')
                old_ax = getappdata(hs.ax(4), 'FreqAxes');
                if ishandle(old_ax), delete(old_ax); end
            end
            
            if show_freq
                % Create two subplots side-by-side WITHIN the 4th quadrant
                % Make plots smaller to show title and labels
                gap = 0.02;
                width = (origPos(3) - 0.08) / 2;
                height = origPos(4) - 0.15; % Leave more room for title and labels
                bottom = origPos(2) + 0.08; % Shift up from bottom
                
                % LEFT: Time domain (real part, no normalization)
                ax_time = axes('Parent', get(hs.ax(4), 'Parent'), 'Units', 'normalized');
                set(ax_time, 'Position', [origPos(1)+0.03 bottom width height]);
                
                % Plot real part of time domain - NO NORMALIZATION
                plot(ax_time, t_ms, real(fid), 'b-', 'LineWidth', 2.0);
                
                % FID-A style formatting with consistent colors
                set(ax_time, 'XColor', 'w', 'YColor', 'w', 'Color', 'w', ...
                    'FontSize', 12, 'LineWidth', 1.5, 'Box', 'off');
                xlabel(ax_time, 'Time (ms)', 'FontSize', 14, 'Color', 'w');
                ylabel(ax_time, '', 'FontSize', 14);
                title(ax_time, sprintf('Time (real) [%d,%d,%d]', ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3)), ...
                    'FontSize', 12, 'Color', 'w');
                xlim(ax_time, [0 max(t_ms)]);
                
                % RIGHT: Frequency domain (real part, no normalization)
                ax_freq = axes('Parent', get(hs.ax(4), 'Parent'), 'Units', 'normalized');
                set(ax_freq, 'Position', [origPos(1)+0.06+width bottom width height]);
                
                % Plot real part of frequency domain - NO NORMALIZATION
                if isfield(vox_MRS, 'ppm') && isfield(vox_MRS, 'specs')
                    ppm = vox_MRS.ppm;
                    specs = real(vox_MRS.specs); % real PART ONLY
                    
                    plot(ax_freq, ppm, specs, 'b-', 'LineWidth', 2.0);
                    
                    % FID-A style formatting with consistent colors
                    set(ax_freq, 'XDir', 'reverse', 'XColor', 'w', 'YColor', 'w', ...
                        'Color', 'w', 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'off');
                    xlabel(ax_freq, 'Chemical Shift (ppm)', 'FontSize', 14, 'Color', 'w');
                    ylabel(ax_freq, '', 'FontSize', 14);
                    title(ax_freq, sprintf('Spectrum (real) [%d,%d,%d]', ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3)), ...
                        'FontSize', 12, 'Color', 'w');
                    xlim(ax_freq, [0 6]); % Typical metabolite range
                    
                    % Add metabolite reference lines
                    hold(ax_freq, 'on');
                    y_lim = ylim(ax_freq);
                    metabolites = struct('NAA', 2.01, 'Cr', 3.03, 'Cho', 3.22);
                    met_names = fieldnames(metabolites);
                    % for i = 1:length(met_names)
                    %     met_ppm = metabolites.(met_names{i});
                    %     if met_ppm >= 0 && met_ppm <= 6
                    %         line(ax_freq, [met_ppm met_ppm], y_lim, ...
                    %             'Color', [0.7 0.7 0.7], 'LineStyle', '--', 'LineWidth', 0.5);
                    %         text(ax_freq, met_ppm, y_lim(2)*0.98, met_names{i}, ...
                    %             'FontSize', 9, 'HorizontalAlignment', 'center', ...
                    %             'VerticalAlignment', 'top', 'Color', [0.5 0.5 0.5]);
                    %     end
                    % end
                    hold(ax_freq, 'off');
                    
                    fprintf('Displaying real part of spectrum (not magnitude)\n');
                else
                    % Fallback if ppm/specs not available
                    text(ax_freq, 0.5, 0.5, 'Error: Missing ppm or specs field', ...
                        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                        'FontSize', 12, 'Color', 'w');
                end
                
                % Store axes handles
                setappdata(hs.ax(4), 'TimeAxes', ax_time);
                setappdata(hs.ax(4), 'FreqAxes', ax_freq);
                
            else
                % Only time domain available - use full width
                cla(hs.ax(4));
                height = origPos(4) - 0.15;
                bottom = origPos(2) + 0.08;
                set(hs.ax(4), 'Position', [origPos(1)+0.03 bottom origPos(3)-0.06 height]);
                
                % Plot real part of time domain - NO NORMALIZATION
                plot(hs.ax(4), t_ms, real(fid), 'b-', 'LineWidth', 2.0);
                
                % FID-A style formatting with consistent colors
                set(hs.ax(4), 'XColor', 'w', 'YColor', 'w', 'Color', 'w', ...
                    'FontSize', 12, 'LineWidth', 1.5, 'Box', 'off');
                xlabel(hs.ax(4), 'Time (ms)', 'FontSize', 14, 'Color', 'w');
                ylabel(hs.ax(4), '', 'FontSize', 14);
                title(hs.ax(4), sprintf('Time Domain (real) [%d,%d,%d]', ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3)), ...
                    'FontSize', 12, 'Color', 'w');
                xlim(hs.ax(4), [0 max(t_ms)]);
                
                fprintf('Displaying real part of FID (not magnitude)\n');
            end
            
            % Hide colorbar when showing spectrum
            if ishandle(hs.colorbar) && isvalid(hs.colorbar)
                set(hs.colorbar, 'Visible', 'off');
            end
            
            fprintf('Spectrum displayed successfully in 4th quadrant\n');
        else
            fprintf('ERROR: ax(4) is not a valid handle\n');
        end
        
    catch ME
        fprintf('ERROR in spectrum display: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
        end
    end
    
    fprintf('=== End spectrum display ===\n');
end
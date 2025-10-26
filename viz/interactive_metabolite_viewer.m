function interactive_metabolite_viewer(ftSpec_data, t1_nifti_path, varargin)
% INTERACTIVE_METABOLITE_VIEWER - WORKING VERSION with visible T1 scaling
%
% T1 is scaled AND padded so the canvas size stays the same
% This makes alignment changes actually visible!
%
% USAGE:
%   interactive_metabolite_viewer(ftSpec_rmlip, 'path/to/T1.nii.gz', 'sliceNumber', 141, 't1ScaleFactor', 0.8)

fprintf('\n========================================\n');
fprintf('METABOLITE VIEWER - T1 SCALING ENABLED\n');
fprintf('========================================\n\n');

%% Parse inputs
p = inputParser;
addRequired(p, 'ftSpec_data');
addRequired(p, 't1_nifti_path');
addParameter(p, 'sliceNumber', [], @isnumeric);
addParameter(p, 't1ScaleFactor', 0.8, @isnumeric);
parse(p, ftSpec_data, t1_nifti_path, varargin{:});
opts = p.Results;

%% Load T1 data
fprintf('[1/4] Loading T1...\n');
try
    t1_data = niftiread(opts.t1_nifti_path);
catch
    V = spm_vol(char(opts.t1_nifti_path));
    t1_data = spm_read_vols(V);
end

%% Determine slice
fprintf('[2/4] Determining slice...\n');
if ~isempty(opts.sliceNumber)
    slice_number = opts.sliceNumber;
else
    slice_number = round(size(t1_data, 3) / 2);
end
fprintf('  Using slice: %d\n', slice_number);

%% Extract MRSI data
fprintf('[3/4] Extracting MRSI data...\n');
mrsi_sz_x = ftSpec_data.sz(ftSpec_data.dims.x);
mrsi_sz_y = ftSpec_data.sz(ftSpec_data.dims.y);

if isfield(ftSpec_data, 'voxelSize')
    voxel_size_x_mm = ftSpec_data.voxelSize.x;
    voxel_size_y_mm = ftSpec_data.voxelSize.y;
else
    voxel_size_x_mm = 10;
    voxel_size_y_mm = 10;
end

ppm = ftSpec_data.ppm;
spec_data_full = real(ftSpec_data.data);
fprintf('  Grid: %d x %d voxels (%.1f x %.1f mm)\n', mrsi_sz_x, mrsi_sz_y, voxel_size_x_mm, voxel_size_y_mm);

%% Prepare T1 slice with scaling and padding
fprintf('[4/4] Preparing T1 with %.0f%% scaling...\n', opts.t1ScaleFactor*100);

t1_slice_orig = rot90(double(t1_data(:, :, slice_number)));
[orig_rows, orig_cols] = size(t1_slice_orig);
fprintf('  Original T1: %d x %d pixels\n', orig_rows, orig_cols);

% Scale T1
scale_factor = opts.t1ScaleFactor;
new_rows = round(orig_rows * scale_factor);
new_cols = round(orig_cols * scale_factor);

[X_orig, Y_orig] = meshgrid(1:orig_cols, 1:orig_rows);
[X_new, Y_new] = meshgrid(linspace(1, orig_cols, new_cols), linspace(1, orig_rows, new_rows));
t1_scaled = interp2(X_orig, Y_orig, t1_slice_orig, X_new, Y_new, 'linear');

fprintf('  Scaled T1: %d x %d pixels (%.0f%%)\n', new_rows, new_cols, scale_factor*100);

% Add padding to maintain original canvas size
pad_rows = (orig_rows - new_rows) / 2;
pad_cols = (orig_cols - new_cols) / 2;

% Create padded image with black background
t1_slice = zeros(orig_rows, orig_cols);
row_start = round(pad_rows) + 1;
row_end = row_start + new_rows - 1;
col_start = round(pad_cols) + 1;
col_end = col_start + new_cols - 1;

t1_slice(row_start:row_end, col_start:col_end) = t1_scaled;

fprintf('  Final T1 (with padding): %d x %d pixels\n', orig_rows, orig_cols);
fprintf('  T1 is centered with black padding around it\n');

% Calculate voxel layout on ORIGINAL canvas size
t1_rows = orig_rows;
t1_cols = orig_cols;
pixels_per_x = t1_cols / mrsi_sz_x;
pixels_per_y = t1_rows / mrsi_sz_y;
voxel_centers_x = ((0:mrsi_sz_x-1) + 0.5) * pixels_per_x;
voxel_centers_y = ((0:mrsi_sz_y-1) + 0.5) * pixels_per_y;

%% Create GUI
main_fig = figure('Position', [50, 50, 1600, 800], 'Color', 'white', ...
                  'Name', sprintf('Metabolite Viewer - T1 Scale: %.0f%%', scale_factor*100), ...
                  'NumberTitle', 'off');

control_panel = uipanel('Parent', main_fig, 'Position', [0.02, 0.80, 0.63, 0.18], ...
                        'Title', 'Controls', 'FontSize', 12, 'FontWeight', 'bold');

uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'Metabolite:', ...
          'FontSize', 11, 'Units', 'normalized', 'Position', [0.02, 0.5, 0.12, 0.35], ...
          'HorizontalAlignment', 'right');

metabolite_dropdown = uicontrol('Parent', control_panel, 'Style', 'popupmenu', ...
          'String', {'Full Spectrum', 'NAA (2.0 ppm)', 'Creatine (3.0 ppm)', 'Choline (3.2 ppm)', ...
                     'Lactate (1.3 ppm)', 'Lipids (0.9-1.4 ppm)', 'Custom Range...'}, ...
          'FontSize', 11, 'Units', 'normalized', 'Position', [0.15, 0.5, 0.20, 0.35]);

uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'PPM Min:', ...
          'FontSize', 11, 'Units', 'normalized', 'Position', [0.37, 0.5, 0.10, 0.35], ...
          'HorizontalAlignment', 'right');
ppm_min_edit = uicontrol('Parent', control_panel, 'Style', 'edit', 'String', '0.2', ...
          'FontSize', 11, 'Units', 'normalized', 'Position', [0.48, 0.5, 0.08, 0.35]);

uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'PPM Max:', ...
          'FontSize', 11, 'Units', 'normalized', 'Position', [0.58, 0.5, 0.10, 0.35], ...
          'HorizontalAlignment', 'right');
ppm_max_edit = uicontrol('Parent', control_panel, 'Style', 'edit', 'String', '4.3', ...
          'FontSize', 11, 'Units', 'normalized', 'Position', [0.69, 0.5, 0.08, 0.35]);

update_button = uicontrol('Parent', control_panel, 'Style', 'pushbutton', 'String', 'Update', ...
          'FontSize', 12, 'FontWeight', 'bold', 'Units', 'normalized', ...
          'Position', [0.80, 0.5, 0.15, 0.35], 'BackgroundColor', [0.2 0.8 0.2]);

uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'Line Width:', ...
          'FontSize', 10, 'Units', 'normalized', 'Position', [0.02, 0.05, 0.12, 0.30], ...
          'HorizontalAlignment', 'right');
linewidth_slider = uicontrol('Parent', control_panel, 'Style', 'slider', ...
          'Min', 0.5, 'Max', 3, 'Value', 1.2, 'Units', 'normalized', ...
          'Position', [0.15, 0.05, 0.15, 0.30]);

uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'Y Scale:', ...
          'FontSize', 10, 'Units', 'normalized', 'Position', [0.32, 0.05, 0.12, 0.30], ...
          'HorizontalAlignment', 'right');
yscale_slider = uicontrol('Parent', control_panel, 'Style', 'slider', ...
          'Min', 0.3, 'Max', 3.0, 'Value', 0.8, 'Units', 'normalized', ...
          'Position', [0.45, 0.05, 0.15, 0.30]);

status_text = uicontrol('Parent', control_panel, 'Style', 'text', ...
          'String', sprintf('T1: %.0f%% scaled | GLOBAL', scale_factor*100), ...
          'FontSize', 10, 'Units', 'normalized', 'Position', [0.62, 0.05, 0.35, 0.30], ...
          'HorizontalAlignment', 'left', 'ForegroundColor', [0 0.5 0]);

overlay_ax = axes('Parent', main_fig, 'Position', [0.05, 0.08, 0.60, 0.68]);
spec_panel = uipanel('Parent', main_fig, 'Position', [0.67, 0.05, 0.31, 0.93], ...
                     'Title', 'Click spectrum', 'FontSize', 12);
spec_ax = axes('Parent', spec_panel, 'Position', [0.12, 0.12, 0.82, 0.80]);

% Store data
data_struct = struct();
data_struct.ftSpec_data = ftSpec_data;
data_struct.ppm = ppm;
data_struct.spec_data_full = spec_data_full;
data_struct.t1_slice = t1_slice;
data_struct.mrsi_sz_x = mrsi_sz_x;
data_struct.mrsi_sz_y = mrsi_sz_y;
data_struct.voxel_size_x_mm = voxel_size_x_mm;
data_struct.pixels_per_x = pixels_per_x;
data_struct.pixels_per_y = pixels_per_y;
data_struct.voxel_centers_x = voxel_centers_x;
data_struct.voxel_centers_y = voxel_centers_y;
data_struct.overlay_ax = overlay_ax;
data_struct.spec_ax = spec_ax;
data_struct.ppm_min_edit = ppm_min_edit;
data_struct.ppm_max_edit = ppm_max_edit;
data_struct.linewidth_slider = linewidth_slider;
data_struct.yscale_slider = yscale_slider;
data_struct.status_text = status_text;
data_struct.metabolite_dropdown = metabolite_dropdown;

guidata(main_fig, data_struct);

set(metabolite_dropdown, 'Callback', {@metabolite_changed, main_fig});
set(update_button, 'Callback', {@update_spectra, main_fig});
set(linewidth_slider, 'Callback', {@update_spectra, main_fig});
set(yscale_slider, 'Callback', {@update_spectra, main_fig});

update_spectra([], [], main_fig);

fprintf('\n========================================\n');
fprintf('READY! Try different scale factors:\n');
fprintf('  0.70 = very small T1\n');
fprintf('  0.80 = default (current)\n');
fprintf('  0.90 = slightly smaller T1\n');
fprintf('========================================\n\n');

end

function metabolite_changed(src, ~, fig_handle)
    data_struct = guidata(fig_handle);
    ranges = {[0.2, 4.3]; [1.9, 2.1]; [2.9, 3.1]; [3.1, 3.3]; [1.2, 1.4]; [0.8, 1.5]; [0.2, 4.3]};
    selection = get(src, 'Value');
    if selection <= length(ranges) && selection ~= 7
        set(data_struct.ppm_min_edit, 'String', num2str(ranges{selection}(1)));
        set(data_struct.ppm_max_edit, 'String', num2str(ranges{selection}(2)));
        update_spectra([], [], fig_handle);
    end
end

function update_spectra(~, ~, fig_handle)
    data_struct = guidata(fig_handle);
    
    ppm_min = str2double(get(data_struct.ppm_min_edit, 'String'));
    ppm_max = str2double(get(data_struct.ppm_max_edit, 'String'));
    linewidth = get(data_struct.linewidth_slider, 'Value');
    yscale = get(data_struct.yscale_slider, 'Value');
    
    if isnan(ppm_min) || isnan(ppm_max) || ppm_min >= ppm_max, return; end
    
    ppm_full = data_struct.ppm;
    spec_data_full = data_struct.spec_data_full;
    t1_slice = data_struct.t1_slice;
    mrsi_sz_x = data_struct.mrsi_sz_x;
    mrsi_sz_y = data_struct.mrsi_sz_y;
    voxel_size_x_mm = data_struct.voxel_size_x_mm;
    pixels_per_x = data_struct.pixels_per_x;
    pixels_per_y = data_struct.pixels_per_y;
    voxel_centers_x = data_struct.voxel_centers_x;
    voxel_centers_y = data_struct.voxel_centers_y;
    overlay_ax = data_struct.overlay_ax;
    
    ppm_mask = (ppm_full >= ppm_min) & (ppm_full <= ppm_max);
    ppm_subset = ppm_full(ppm_mask);
    spec_data = spec_data_full(ppm_mask, :, :);
    
    % GLOBAL SCALING
    global_spectrum_height = max(max(spec_data, [], 1) - min(spec_data, [], 1), [], 'all');
    if global_spectrum_height > 0
        global_scale_factor = yscale * (0.8 * pixels_per_y) / global_spectrum_height;
    else
        global_scale_factor = 1;
    end
    
    cla(overlay_ax);
    axes(overlay_ax);
    imagesc(overlay_ax, t1_slice);
    colormap(overlay_ax, 'gray');
    axis(overlay_ax, 'image', 'off');
    hold(overlay_ax, 'on');
    
    blue_color = [0, 0, 1];
    
    n_plotted = 0;
    for ix = 1:mrsi_sz_x
        for iy = 1:mrsi_sz_y
            spectrum = spec_data(:, iy, ix);
            if all(spectrum == 0) || all(isnan(spectrum)), continue; end
            
            spectrum_scaled = spectrum * global_scale_factor;
            spectrum_width_fraction = min(0.9, 0.5 + (voxel_size_x_mm / 20));
            
            if length(ppm_subset) > 1
                ppm_norm = (ppm_subset - min(ppm_subset)) / (max(ppm_subset) - min(ppm_subset));
            else
                ppm_norm = 0.5;
            end
            
            x_coords = voxel_centers_x(ix) + (ppm_norm - 0.5) * pixels_per_x * spectrum_width_fraction;
            y_coords = voxel_centers_y(iy) + spectrum_scaled;
            
            plot(overlay_ax, x_coords, y_coords, 'Color', blue_color, 'LineWidth', linewidth, ...
                 'ButtonDownFcn', {@show_spectrum, fig_handle, ix, iy});
            n_plotted = n_plotted + 1;
        end
    end
    
    for i = 1:mrsi_sz_x-1
        plot(overlay_ax, [i*pixels_per_x, i*pixels_per_x], [1, size(t1_slice,1)], 'y:', 'LineWidth', 0.5);
    end
    for j = 1:mrsi_sz_y-1
        plot(overlay_ax, [1, size(t1_slice,2)], [j*pixels_per_y, j*pixels_per_y], 'y:', 'LineWidth', 0.5);
    end
    
    title(overlay_ax, sprintf('PPM: %.2f-%.2f | %d spectra | GLOBAL', ppm_min, ppm_max, n_plotted), ...
          'FontSize', 12, 'FontWeight', 'bold');
    hold(overlay_ax, 'off');
    
    set(data_struct.status_text, 'String', sprintf('%d spectra | GLOBAL', n_plotted), ...
        'ForegroundColor', [0 0.5 0]);
end

function show_spectrum(~, ~, fig_handle, voxel_x, voxel_y)
    data_struct = guidata(fig_handle);
    full_ppm = data_struct.ftSpec_data.ppm;
    full_spec = real(data_struct.ftSpec_data.data(:, voxel_y, voxel_x));
    spec_ax = data_struct.spec_ax;
    
    cla(spec_ax);
    plot(spec_ax, full_ppm, full_spec, 'b-', 'LineWidth', 2);
    set(spec_ax, 'XDir', 'reverse');
    grid(spec_ax, 'on');
    xlabel(spec_ax, 'Chemical Shift (ppm)', 'FontSize', 11);
    ylabel(spec_ax, 'Signal', 'FontSize', 11);
    title(spec_ax, sprintf('Voxel [%d, %d]', voxel_x, voxel_y), 'FontSize', 12, 'FontWeight', 'bold');
end
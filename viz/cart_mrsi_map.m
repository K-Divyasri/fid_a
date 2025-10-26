function out = convert_ftSpec_with_spec2nii_exact_positioning(S_time, S_freq, outPath, refT1Path, emptyNiiPath, opts)
% Convert using spec2nii CENTER and Z orientation, YOUR data for XY voxel sizes

if nargin < 6 || isempty(opts), opts = struct; end
def = struct('ZOffsetMM', 0, 'ManualOffset', [0,0,0], 'DebugMode', true, ...
             'MakeFirstPointOverlay', true, 'OverlayIndex', 1, ...
             'OverlayOutPath', '');
fn = fieldnames(def);
for k=1:numel(fn), if ~isfield(opts,fn{k}), opts.(fn{k}) = def.(fn{k}); end, end

assert(exist('nii_tool','file')==2, 'Add dicm2nii to path (for nii_tool).');

% Validate time domain data (ccav_w)
assert(isfield(S_time,'dims') && S_time.dims.x>0 && S_time.dims.y>0, 'S_time.dims.x/y must be set.');
assert(isfield(S_time,'data') && ndims(S_time.data)==3, 'S_time.data must be [T x X x Y].');
assert(isfield(S_time,'affineMatrix') && all(size(S_time.affineMatrix)==[4,4]), 'S_time.affineMatrix 4x4 required.');

% Validate frequency domain data
assert(isfield(S_freq,'dims') && S_freq.dims.x>0 && S_freq.dims.y>0, 'S_freq.dims.x/y must be set.');
assert(isfield(S_freq,'data') && ndims(S_freq.data)==3, 'S_freq.data must be [F x X x Y].');

% Time dimension index
if isfield(S_time.dims,'t') && S_time.dims.t>0
    dSpec_time = S_time.dims.t;
elseif isfield(S_time.dims,'f') && S_time.dims.f>0
    dSpec_time = S_time.dims.f;
else
    error('S_time needs dims.t or dims.f > 0');
end

% ---- Reformat TIME domain ----
permOrder_time = [S_time.dims.x S_time.dims.y dSpec_time];
vol_time = permute(S_time.data, permOrder_time);
[X, Y, T] = size(vol_time);
vol_time = reshape(vol_time, [X Y 1 T]);
vol_time = single(vol_time);

fprintf('\n=== MRSI CONVERSION ===\n');
fprintf('Grid: %d x %d x 1 x %d points\n', X, Y, T);

% ---- Extract XY voxel sizes from YOUR data ----
dx_data = norm(S_time.affineMatrix(1:3,1));
dy_data = norm(S_time.affineMatrix(1:3,2));

fprintf('YOUR data XY voxel sizes: [%.2f, %.2f] mm\n', dx_data, dy_data);

% ---- Load spec2nii ----
fprintf('\n=== LOADING SPEC2NII ===\n');
fprintf('Path: %s\n', emptyNiiPath);
assert(exist(emptyNiiPath,'file')==2, 'spec2nii not found');

Ht = nii_tool('hdr', emptyNiiPath);
fprintf('spec2nii dimensions: [%d %d %d %d]\n', Ht.dim(2:5));

assert(Ht.sform_code>0, 'spec2nii has no sform');
if isfield(Ht, 'sform_mat') && ~isempty(Ht.sform_mat)
    sform_spec2nii = [double(Ht.sform_mat); 0 0 0 1];
else
    sform_spec2nii = [double(Ht.srow_x); double(Ht.srow_y); double(Ht.srow_z); 0 0 0 1];
end

fprintf('spec2nii sform:\n'); disp(sform_spec2nii);

% ---- Extract orientation vectors from spec2nii ----
xd_spec2nii = sform_spec2nii(1:3,1);
yd_spec2nii = sform_spec2nii(1:3,2);
zd_spec2nii = sform_spec2nii(1:3,3);

x_dir_spec2nii = xd_spec2nii / norm(xd_spec2nii);
y_dir_spec2nii = yd_spec2nii / norm(yd_spec2nii);
z_dir_spec2nii = zd_spec2nii / norm(zd_spec2nii);

dx_spec2nii = norm(xd_spec2nii);
dy_spec2nii = norm(yd_spec2nii);
dz_spec2nii = norm(zd_spec2nii);

fprintf('spec2nii orientations:\n');
fprintf('  X_dir: [%.4f, %.4f, %.4f], size: %.2f mm\n', x_dir_spec2nii, dx_spec2nii);
fprintf('  Y_dir: [%.4f, %.4f, %.4f], size: %.2f mm\n', y_dir_spec2nii, dy_spec2nii);
fprintf('  Z_dir: [%.4f, %.4f, %.4f], size: %.2f mm\n', z_dir_spec2nii, dz_spec2nii);

% ---- Extract CENTER from spec2nii ----
dims_spec2nii = double(Ht.dim(2:4));
center_vox_spec2nii = [(dims_spec2nii(1)-1)/2; (dims_spec2nii(2)-1)/2; 0; 1];
center_mm_spec2nii = sform_spec2nii * center_vox_spec2nii;
spec2nii_center = center_mm_spec2nii(1:3);

fprintf('spec2nii center: [%.2f, %.2f, %.2f] mm\n', spec2nii_center);

% ---- FINAL: spec2nii orientation + center, YOUR XY voxel sizes ----
x_dir_final = x_dir_spec2nii;  % spec2nii X orientation
y_dir_final = y_dir_spec2nii;  % spec2nii Y orientation
z_dir_final = z_dir_spec2nii;  % spec2nii Z orientation

dx_final = dx_data;  % YOUR X voxel size
dy_final = dy_data;  % YOUR Y voxel size
dz_final = dz_spec2nii;  % spec2nii Z thickness

fprintf('\n=== FINAL CONFIG ===\n');
fprintf('Orientations: FROM SPEC2NII\n');
fprintf('  X_dir: [%.4f, %.4f, %.4f], size: %.2f mm (YOUR size)\n', x_dir_final, dx_final);
fprintf('  Y_dir: [%.4f, %.4f, %.4f], size: %.2f mm (YOUR size)\n', y_dir_final, dy_final);
fprintf('  Z_dir: [%.4f, %.4f, %.4f], size: %.2f mm (spec2nii)\n', z_dir_final, dz_final);
fprintf('Center: FROM SPEC2NII [%.2f, %.2f, %.2f] mm\n', spec2nii_center);

% ---- Apply offsets to spec2nii center ----
target_center = spec2nii_center + [opts.ManualOffset(1); opts.ManualOffset(2); opts.ManualOffset(3) + opts.ZOffsetMM];

if norm([opts.ManualOffset(:); opts.ZOffsetMM]) > 0
    fprintf('Applied offsets: Manual [%.2f, %.2f, %.2f] + Z %.2f mm\n', opts.ManualOffset, opts.ZOffsetMM);
    fprintf('Final center: [%.2f, %.2f, %.2f] mm\n', target_center);
end

% ---- Build affine ----
x_vec = x_dir_final * dx_final;
y_vec = y_dir_final * dy_final;
z_vec = z_dir_final * dz_final;

halfX = ((X-1)/2) * x_vec;
halfY = ((Y-1)/2) * y_vec;
corner_pos = target_center - halfX - halfY;

A_final = eye(4);
A_final(1:3,1) = x_vec;
A_final(1:3,2) = y_vec;
A_final(1:3,3) = z_vec;
A_final(1:3,4) = corner_pos;

if opts.DebugMode
    fprintf('\n=== AFFINE MATRIX ===\n');
    disp(A_final);
    
    % Verify center
    center_vox = [(X-1)/2; (Y-1)/2; 0; 1];
    computed = A_final * center_vox;
    fprintf('\nVerification:\n');
    fprintf('  Target:   [%.2f, %.2f, %.2f] mm\n', target_center);
    fprintf('  Computed: [%.2f, %.2f, %.2f] mm\n', computed(1:3));
    fprintf('  Error:    %.4f mm\n', norm(computed(1:3)-target_center));
    
    % Grid coverage
    fprintf('\nGrid coverage:\n');
    fprintf('  X: %.2f mm (%d voxels × %.2f mm)\n', X*dx_final, X, dx_final);
    fprintf('  Y: %.2f mm (%d voxels × %.2f mm)\n', Y*dy_final, Y, dy_final);
end

% ---- T1 distance ----
if ~isempty(refT1Path) && exist(refT1Path,'file')
    T1c = local_get_center(refT1Path);
    if ~isempty(T1c)
        fprintf('\nT1 center: [%.2f, %.2f, %.2f] mm\n', T1c);
        fprintf('Distance from T1: %.2f mm\n', norm(target_center-T1c));
    end
end

% ---- Write NIfTI ----
fprintf('\n=== WRITING NIFTI ===\n');
nii_time = nii_tool('init', vol_time);
nii_time.hdr.dim(1)      = int16(4);
nii_time.hdr.dim(2:5)    = int16([X Y 1 T]);
nii_time.hdr.pixdim(2:5) = single([dx_final dy_final dz_final 5e-6]);
nii_time.hdr.xyzt_units  = uint8(2+8);
nii_time.hdr.datatype    = int16(32);
nii_time.hdr.bitpix      = int16(64);
nii_time.hdr.sform_code  = int16(1);

if isfield(nii_time.hdr,'sform_mat')
    nii_time.hdr.sform_mat = single(A_final(1:3,:));
else
    nii_time.hdr.srow_x = single(A_final(1,:));
    nii_time.hdr.srow_y = single(A_final(2,:));
    nii_time.hdr.srow_z = single(A_final(3,:));
end
nii_time.hdr.qform_code  = int16(0);
nii_time.hdr.descrip     = sprintf('MRSI [%dx%dx1x%d], spec2nii pose', X, Y, T);
nii_time.hdr.intent_code = int16(2006);
nii_time.hdr.intent_name = 'MRSI_TIME';

nii_tool('save', nii_time, outPath);
fprintf('Saved: %s\n', outPath);

% ---- Save ftSpec ----
[p, n, ~] = fileparts(outPath);
if endsWith(n, '.nii'), n = extractBefore(n, '.nii'); end
ftSpecPath = fullfile(p, [n '_ftSpec.mat']);
ftSpec_smooth_w = S_freq;
save(ftSpecPath, 'ftSpec_smooth_w', '-v7.3');
fprintf('Saved: %s\n', ftSpecPath);

out = struct('mrsi4D_time', outPath, 'ftSpec_mat', ftSpecPath, 'overlay3D', '');

% ---- Overlay ----
if opts.MakeFirstPointOverlay
    idx = max(1, min(T, round(opts.OverlayIndex)));
    mag2d = abs(vol_time(:,:,1,idx));
    mmax = max(mag2d(:));
    if mmax>0, mag2d = mag2d / mmax; end
    mag3d = reshape(mag2d, [X Y 1]);

    niiO = nii_tool('init', single(mag3d));
    niiO.hdr.dim(1)      = int16(3);
    niiO.hdr.dim(2:5)    = int16([X Y 1 1]);
    niiO.hdr.pixdim(2:5) = single([dx_final dy_final dz_final 1]);
    niiO.hdr.xyzt_units  = uint8(2+8);
    niiO.hdr.datatype    = int16(16);
    niiO.hdr.bitpix      = int16(32);
    niiO.hdr.sform_code  = int16(1);
    
    if isfield(niiO.hdr,'sform_mat')
        niiO.hdr.sform_mat = single(A_final(1:3,:));
    else
        niiO.hdr.srow_x = single(A_final(1,:));
        niiO.hdr.srow_y = single(A_final(2,:));
        niiO.hdr.srow_z = single(A_final(3,:));
    end
    niiO.hdr.qform_code = int16(0);

    if isempty(opts.OverlayOutPath)
        overlayPath = fullfile(p, [n '_firstpt_overlay.nii.gz']);
    else
        overlayPath = opts.OverlayOutPath;
    end
    
    nii_tool('save', niiO, overlayPath);
    fprintf('Saved overlay: %s\n', overlayPath);
    out.overlay3D = overlayPath;
end

fprintf('==============================================\n');
fprintf('DONE: spec2nii center + orientations, YOUR XY sizes\n');
fprintf('==============================================\n');
end

function c = local_get_center(niiPath)
c = [];
try
    H = nii_tool('hdr', niiPath);
    if H.sform_code>0
        sz = double(H.dim(2:4));
        vox = (sz-1)/2;
        if isfield(H,'sform_mat') && ~isempty(H.sform_mat)
            A = [double(H.sform_mat); 0 0 0 1];
        else
            A = [double(H.srow_x); double(H.srow_y); double(H.srow_z); 0 0 0 1];
        end
        w = A * [vox'; 1];
        c = w(1:3);
    end
catch
end
end
% After your rosette pipeline:
t1Fixed = 'C:\Users\divya\Downloads\Lexar\DataForDivya\invivo_10Sep2024_40x40_rose\invivo_10Sep2024_40x40_rose\output\T1_with_sform.nii.gz';
mrsiOut = 'mrsi_time_domain.nii.gz';
emptyNiiPath = "C:\Users\divya\Downloads\phantom_shift\empty niftis\meas_MID00520_FID124507_head_csi_fid_24x24_te2p3_tr750_w_empty.nii.gz";

% Convert: ccav_w (time, for display) + ftSpec_smooth_w (freq, for FID-A)
out = convert_ftSpec_with_spec2nii_exact_positioning(...
    ccav_w, ...              % Time domain (left plot)
    ftSpec_smooth, ...     % Frequency domain (for op_CSItoMRS + op_plotspec)
    mrsiOut, ...
    t1Fixed, ...
    emptyNiiPath, ...
    struct('ZOffsetMM', 0, 'ManualOffset', [0 0 0], ...
           'MakeFirstPointOverlay', true));

% View - will auto-load ftSpec_smooth_w and use FID-A functions
nii_viewer(t1Fixed, out.mrsi4D_time);
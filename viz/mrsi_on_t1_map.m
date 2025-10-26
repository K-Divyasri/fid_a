function out = mrsi_on_t1_map(S_time, S_freq, outPath, refT1Path, emptyNiiPath, opts)
% Convert ccav_w (time) and ftSpec_smooth_w (frequency) using EXACT positioning
% S_time = ccav_w (time domain for display)
% S_freq = ftSpec_smooth_w (frequency domain for FID-A processing)

if nargin < 6 || isempty(opts), opts = struct; end
def = struct('ZOffsetMM', 0, 'ManualOffset', [0,0,0], 'DebugMode', true, ...
             'MakeFirstPointOverlay', true, 'OverlayIndex', 1, ...
             'OverlayOutPath', '', 'FreqOutPath', '');
fn = fieldnames(def);
for k=1:numel(fn), if ~isfield(opts,fn{k}), opts.(fn{k}) = def.(fn{k}); end, end

assert(exist('nii_tool','file')==2, 'Add dicm2nii to path (for nii_tool).');

% Validate time domain data (ccav_w)
assert(isfield(S_time,'dims') && S_time.dims.x>0 && S_time.dims.y>0, 'S_time.dims.x/y must be set.');
assert(isfield(S_time,'data') && ndims(S_time.data)==3, 'S_time.data must be [T x X x Y].');
assert(isfield(S_time,'affineMatrix') && all(size(S_time.affineMatrix)==[4,4]), 'S_time.affineMatrix 4x4 required.');

% Validate frequency domain data (ftSpec_smooth_w)
assert(isfield(S_freq,'dims') && S_freq.dims.x>0 && S_freq.dims.y>0, 'S_freq.dims.x/y must be set.');
assert(isfield(S_freq,'data') && ndims(S_freq.data)==3, 'S_freq.data must be [F x X x Y].');

% Time dimension index in S_time.data
if isfield(S_time.dims,'t') && S_time.dims.t>0
    dSpec_time = S_time.dims.t;
elseif isfield(S_time.dims,'f') && S_time.dims.f>0
    dSpec_time = S_time.dims.f;
else
    error('S_time needs dims.t or dims.f > 0');
end

% ---- Reformat TIME domain to [X Y 1 T] for NIfTI ----
permOrder_time = [S_time.dims.x S_time.dims.y dSpec_time];
vol_time = permute(S_time.data, permOrder_time);
[X, Y, T] = size(vol_time);
vol_time = reshape(vol_time, [X Y 1 T]);
vol_time = single(vol_time);

fprintf('\n=== MRSI TIME-DOMAIN CONVERSION ===\n');
fprintf('Time domain: %d x %d x 1 x %d points\n', X, Y, T);

% ---- Extract voxel sizes from S_time.affineMatrix ----
dx_ftspec = norm(S_time.affineMatrix(1:3,1));
dy_ftspec = norm(S_time.affineMatrix(1:3,2));
dz_ftspec = norm(S_time.affineMatrix(1:3,3));
fprintf('Voxel sizes: [%.2f, %.2f, %.2f] mm\n', dx_ftspec, dy_ftspec, dz_ftspec);

% ---- Load spec2nii "empty" NIfTI (for pose) ----
fprintf('\n=== LOADING SPEC2NII EXACT POSITIONING ===\n');
fprintf('Empty NIfTI path: %s\n', emptyNiiPath);
assert(exist(emptyNiiPath,'file')==2, 'spec2nii empty NIfTI not found.');

Ht = nii_tool('hdr', emptyNiiPath);
fprintf('spec2nii loaded successfully\n');

assert(Ht.sform_code>0, 'spec2nii NIfTI has no spatial transform (sform_code = 0)');
if isfield(Ht, 'sform_mat') && ~isempty(Ht.sform_mat)
    sform_spec2nii = [double(Ht.sform_mat); 0 0 0 1];
else
    assert(all(isfield(Ht, {'srow_x','srow_y','srow_z'})), 'No sform_mat / srow_* present.');
    sform_spec2nii = [double(Ht.srow_x); double(Ht.srow_y); double(Ht.srow_z); 0 0 0 1];
end

% ---- Center from spec2nii pose ----
dims_spec2nii = double(Ht.dim(2:4));
center_vox = (dims_spec2nii - 1) / 2;
center_mm  = sform_spec2nii * [center_vox'; 1];
slice_center = center_mm(1:3);

% ---- Normalized axes from spec2nii ----
xd = sform_spec2nii(1:3,1); yd = sform_spec2nii(1:3,2); zd = sform_spec2nii(1:3,3);
x_dir = xd / norm(xd);  y_dir = yd / norm(yd);  z_dir = zd / norm(zd);

% ---- Apply offsets ----
target_center = slice_center + [opts.ManualOffset(1); opts.ManualOffset(2); opts.ManualOffset(3) + opts.ZOffsetMM];

fprintf('\n=== TARGET POSITIONING ===\n');
fprintf('Final target center: [%.2f, %.2f, %.2f] mm\n', target_center);

% ---- Scale axes with our voxel sizes ----
x_vec = x_dir * dx_ftspec;
y_vec = y_dir * dy_ftspec;
z_vec = z_dir * dz_ftspec;

halfX = (X-1)/2 * x_vec;
halfY = (Y-1)/2 * y_vec;
corner_pos = target_center - halfX - halfY;

A_final = eye(4);
A_final(1:3,1) = x_vec;
A_final(1:3,2) = y_vec;
A_final(1:3,3) = z_vec;
A_final(1:3,4) = corner_pos;

if opts.DebugMode
    fprintf('\n=== FINAL AFFINE ===\n'); disp(A_final);
end

% ---- Optional: distance to T1 center ----
if ~isempty(refT1Path) && exist(refT1Path,'file')
    T1c = local_get_center(refT1Path);
    if ~isempty(T1c)
        fprintf('Distance from T1 center: %.2f mm\n', norm(target_center-T1c));
    end
end

% ---- Write TIME DOMAIN 4D NIfTI ----
fprintf('\n=== WRITING TIME DOMAIN 4D NIFTI ===\n');
nii_time = nii_tool('init', vol_time);
nii_time.hdr.dim(1)      = int16(4);
nii_time.hdr.dim(2:5)    = int16([X Y 1 T]);
nii_time.hdr.pixdim(2:5) = single([dx_ftspec dy_ftspec dz_ftspec 5e-6]); % Default dwell time
nii_time.hdr.xyzt_units  = uint8(2+8);
nii_time.hdr.datatype    = int16(32); % complex64
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
nii_time.hdr.descrip     = sprintf('MRSI time domain (ccav_w) [%dx%dx1x%d]', X, Y, T);
nii_time.hdr.intent_code = int16(2006);
nii_time.hdr.intent_name = 'MRSI_TIME';

nii_tool('save', nii_time, outPath);
fprintf('Saved time domain: %s\n', outPath);

% ---- Save ftSpec_smooth_w as companion .mat file ----
fprintf('\n=== SAVING ftSpec_smooth_w for FID-A processing ===\n');
[p, n, ~] = fileparts(outPath);
if endsWith(n, '.nii'), n = extractBefore(n, '.nii'); end
ftSpecPath = fullfile(p, [n '_ftSpec.mat']);

ftSpec_smooth_w = S_freq; % Store the full FID-A structure
save(ftSpecPath, 'ftSpec_smooth_w', '-v7.3');
fprintf('Saved ftSpec_smooth_w: %s\n', ftSpecPath);

out = struct('mrsi4D_time', outPath, 'ftSpec_mat', ftSpecPath, 'overlay3D', '');

% ---- Optional: write first-point magnitude overlay ----
if opts.MakeFirstPointOverlay
    idx = max(1, min(T, round(opts.OverlayIndex)));
    mag2d = abs(vol_time(:,:,1,idx));
    
    mmax = max(mag2d(:));
    if mmax>0, mag2d = mag2d / mmax; end
    mag3d = reshape(mag2d, [X Y 1]);

    niiO = nii_tool('init', single(mag3d));
    niiO.hdr.dim(1)      = int16(3);
    niiO.hdr.dim(2:5)    = int16([X Y 1 1]);
    niiO.hdr.pixdim(2:5) = single([dx_ftspec dy_ftspec dz_ftspec 1]);
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
    niiO.hdr.qform_code  = int16(0);
    niiO.hdr.descrip     = sprintf('First-point magnitude overlay (idx=%d)', idx);

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
fprintf('SUMMARY:\n');
fprintf('  Time domain (ccav_w): %s\n', out.mrsi4D_time);
fprintf('  Freq struct (ftSpec): %s\n', out.ftSpec_mat);
if opts.MakeFirstPointOverlay
    fprintf('  Overlay:              %s\n', out.overlay3D);
end
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
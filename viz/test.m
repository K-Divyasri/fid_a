
t1Fixed = 'C:\Users\divya\Downloads\linux-send\fida3 - Copy (2)\data\T1_with_sform.nii.gz';
emptyNiiPath = 'C:\Users\divya\Downloads\linux-send\fida3 - Copy (2)\data\meas_MID01416_FID95677_lb_RosetteSpinEcho_ek_v1a_w_empty.nii.gz';
mrsiOut = 'C:\Users\divya\Downloads\linux-send\fida3 - Copy (2)\data\mrsi_time_domain.nii.gz';

out = mrsi_on_t1_map(...
    ccav_w, ...              % Time domain (left plot)
    ftSpec, ...     % Frequency domain (for op_CSItoMRS + op_plotspec)
    mrsiOut, ...
    t1Fixed, ...
    emptyNiiPath, ...
    struct('ZOffsetMM', 0, 'ManualOffset', [0 0 0], ...
           'MakeFirstPointOverlay', true));

% View - will auto-load ftSpec_smooth_w and use FID-A functions
nii_viewer(t1Fixed, out.mrsi4D_time);

load('C:\Users\divya\Downloads\linux-send\fida3 - Copy (2)\data\map.mat')
load('C:\Users\divya\Downloads\linux-send\fida3 - Copy (2)\data\crlb.mat')
load('C:\Users\divya\Downloads\linux-send\fida3 - Copy (2)\data\LW.mat')
load('C:\Users\divya\Downloads\linux-send\fida3 - Copy (2)\data\SNR.mat')

t1Fixed = 'C:\Users\divya\Downloads\linux-send\fida3 - Copy (2)\data\T1_with_sform.nii.gz';
emptyNiiPath = 'C:\Users\divya\Downloads\linux-send\fida3 - Copy (2)\data\meas_MID01416_FID95677_lb_RosetteSpinEcho_ek_v1a_w_empty.nii.gz';
mrsiOut = 'C:\Users\divya\Downloads\linux-send\fida3 - Copy (2)\data\mrsi_time_domain.nii.gz';
spec2nii_empty_path=emptyNiiPath;
% create_spec2nii_aligned_overlay_viewer(ftSpec_smooth, map, crlb, LW, SNR, t1Fixed, spec2nii_empty_path);

interactive_metabolite_viewer(ftSpec_rmlip, t1Fixed, ...
    'sliceNumber', 141, 't1ScaleFactor', 0.75);
create_spec2nii_aligned_viewer(ftSpec_smooth, map, crlb, LW, SNR, t1Fixed, spec2nii_empty_path, -10);
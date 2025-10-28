
fileName = 'C:\Users\divya\Downloads\linux-send\fid_a - Copy\data\meas_MID00519_FID124506_head_csi_fid_24x24_te2p3_tr750.dat';
fileName_w = 'C:\Users\divya\Downloads\linux-send\fid_a - Copy\data\meas_MID00520_FID124507_head_csi_fid_24x24_te2p3_tr750_w.dat';
kFile ='';
[timeCombined_sh, timeCombined_sh_w] = load_twix2(fileName, fileName_w,kFile);

% K-space to spatial FT
ftSpatial = op_CSIFourierTransform(timeCombined_sh, "", 'spatial', true, 'spectral', false);
ftSpatial_w = op_CSIFourierTransform(timeCombined_sh_w, "", 'spatial', true, 'spectral', false);

% Coil combination (WS first)
[coilCombined_w, phase, weights] = op_CSICombineCoils(ftSpatial_w);
coilCombined = op_CSICombineCoils(ftSpatial, 1, phase, weights);

% Averagetime
ccav = op_CSIAverage(coilCombined);
ccav_w = op_CSIAverage(coilCombined_w);

% Time → Spectral FT
ftSpec = op_CSIFourierTransform(ccav);
ftSpec_w = op_CSIFourierTransform(ccav_w);

% Remove lipids from WS data
ftSpec_rmlip = op_CSIssp(ftSpec, 1.0, 1.88);
ftSpec_rmw = op_CSIRemoveLipids(ftSpec_rmlip, lipidPPMRange=[4.4 5.0], linewidthRange=[0.5 20]);

% [ftSpec_B0corr,ftSpec_B0corr_w,freqMap,R2Map] = op_CSIB0Correction_v2(ftSpec_rmw,ftSpec_w);

% B0 correction using WS as reference\
[ftSpec_B0corr_w, phaseMap, freqMap] = op_CSIB0Correction(ftSpec_w);
ftSpec_B0corr = op_CSIB0Correction(ftSpec_rmw, phaseMap, freqMap);

% Apodization (Gaussian smoothing)
ftSpec_smooth = op_CSIApodize(ftSpec_B0corr, 'functionType', 'gaussian', 'fullWidthHalfMax', 25);
ftSpec_smooth_w = op_CSIApodize(ftSpec_B0corr_w, 'functionType', 'gaussian', 'fullWidthHalfMax', 25);

% Optional plot
op_CSIPlot(ftSpec_smooth);
fprintf('Generating brain and lipid masks...\n');
[brain_area_raw, lipid_ring_raw] = createBrainArea(ftSpec_smooth, ...
    ftSpec_smooth.sz(ftSpec_smooth.dims.x)/2, ...
    ftSpec_smooth.sz(ftSpec_smooth.dims.y)/2);
ftSpec_smooth.mask.brainmasks = brain_area_raw;
ftSpec_smooth.mask.lipmasks = repmat(lipid_ring_raw, [1, 1, 5]);

lcmDir = cd;

fprintf('Writing LCModel files...\n');
io_CSIwritelcm(ftSpec_smooth_w, fullfile(lcmDir, 'lcm_ftSpec_smooth_w'));
io_CSIwritelcm(ftSpec_smooth, fullfile(lcmDir, 'lcm_ftSpec_smooth'));

fprintf('LCModel files saved to: %s\n', lcmDir);
fprintf('*** NOW RUN LCMODEL ON THESE FILES ***\n');
fprintf('After LCModel completes, continue to the next section...\n');

%% === PLOT SPECTRA ===
op_CSIPlot(ftSpec_smooth);
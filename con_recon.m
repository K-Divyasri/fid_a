%% === INPUTS ===
%% ---------------- USER PATHS ----------------
kFile      = 'C:\Users\divya\Downloads\linux-send\fida3\new\concentric_trajectory_kspace_trajectory_circle.txt';
fileName   ='C:\Users\divya\Downloads\jn_concRings_SE\meas_MID01094_FID11147_jn_concRings_SE.dat';
fileName_w ='C:\Users\divya\Downloads\jn_concRings_SE\meas_MID01095_FID11148_jn_concRings_SE_w.dat';

fprintf('=== TWO FILE NON-CARTESIAN MODE ===\n');
cc_w = io_MRSI_load_twix(fileName_w,kFile);
cc = io_MRSI_load_twix(fileName,kFile);

fprintf('  ADC dwell time set to: %.3e s\n\n', cc.adcDwellTime);

%% ---------------- STEP 3: Combine ADC blocks ----------------
fprintf('Step 3: Combining ADC blocks...\n');
timeCombined   = op_CSICombineTime(cc,   'extras');
timeCombined_w = op_CSICombineTime(cc_w, 'extras');
fprintf('  Data size after combine: [%s]\n', num2str(timeCombined.sz));

% % nComp = 8;  % choose, e.g., 8 virtual channels out of your 32
% [timeCombined,  W] = op_CSICoilCompression_SVD(timeCombined_sh,  nComp);
% [timeCombined_w, ~] = op_CSICoilCompression_SVD(timeCombined_sh_w, nComp, false);
% timeCombined_w.coilCompression.W = W;    % overwrite voxel weights
% %  % overwrite voxel weights
%%=======
% 
 % overwrite voxel weights
%%>>>>>>> 4043b67b13dbd02a105840aa769ff77a65d952b4

%% === RESHAPE ===
timeCombined_rs = reshape_twix_data(timeCombined,kFile);
timeCombined_rs_w = reshape_twix_data(timeCombined_w,kFile);

%% === K-SPACE COtRRECTION AND NUFFT ===
fprintf('Applying k-space correction and NUFFT...\n');
dComp = op_CSIPSFCorrection_jn2(timeCombined_rs, kFile, ...
    'dcfMethod','PipeMenon', 'pmIters',25, 'pmKernel','gaussian', ...
    'modelType','Gaussian', 'sigma',0.45, 'isPlotWeights',true, 'isPlotRosette',false);
dComp_w = op_CSIPSFCorrection_jn2(timeCombined_rs_w, kFile, ...
    'dcfMethod','PipeMenon', 'pmIters',25, 'pmKernel','gaussian', ...
    'modelType','Gaussian', 'sigma',0.45, 'isPlotWeights',true, 'isPlotRosette',false);


ftSpatial = op_CSIFourierTransform(dComp, kFile, 'spatial', true, 'spectral', false);
ftSpatial_w = op_CSIFourierTransform(dComp_w, kFile, 'spatial', true, 'spectral', false);

%% === COIL COMBINATION ===
fprintf('Combining coils...\n');
[coilCombined_w, phase, weights] = op_CSICombineCoils1(ftSpatial_w);
coilCombined = op_CSICombineCoils1(ftSpatial, 1, phase, weights);

%% === COMBINE AVERAGES ===
fprintf('Combining averages...\n');
ccav = op_CSIAverage(coilCombined);
ccav_w = op_CSIAverage(coilCombined_w);

%% === SPECTRAL FT ===
fprintf('Performing spectral Fourier transform...\n');
ftSpec = op_CSIFourierTransform(ccav);
ftSpec_w = op_CSIFourierTransform(ccav_w);

% ftSpec.dims.t=1;
% ftSpec.dims.f=0;
% ftSpec_w.dims.t=1;
% ftSpec_w.dims.f=0;
%Remove the residual lipids from the water suppressed data:
ftSpec_rmlip = op_CSIssp(ftSpec,0.8,1.88);

%Remove  residual water from the water suppressed data:
ftSpec_rmw = op_CSIRemoveLipids(ftSpec_rmlip,lipidPPMRange=[4.4 5.0],linewidthRange=[0.5 20]);

%Do a B0 correction:
[ftSpec_B0corr_w,phaseMap,freqMap] = op_CSIB0Correction(ftSpec_w);
[ftSpec_B0corr] = op_CSIB0Correction(ftSpec_rmw,phaseMap,freqMap);

%BEFORE Spatial smoothing, we should try to apply a brain mask so that
%the voxels with mostly lipids are removed.  This way they will not
%bleed into the brain following spatial smoothing.   But this is for
%another day.  

%Finally, do spatial smoothing
[ftSpec_smooth] = op_CSIApodize(ftSpec_B0corr,'functionType','gaussian','fullWidthHalfMax',15);
[ftSpec_smooth_w] = op_CSIApodize(ftSpec_B0corr_w,'functionType','gaussian','fullWidthHalfMax',15);

% %% === LIPID REMOVAL ===
% fprintf('Removing lipids...\n');
% ftSpec_rmlip = op_CSIssp(ftSpec, 1.0, 1.88);
% ftSpec_rmw = op_CSIRemoveLipids(ftSpec_rmlip, lipidPPMRange=[4.4 5.0], linewidthRange=[0.5 20]);
% 
% 
% %Do a B0 correction:
% [ftSpec_B0corr_w,phaseMap,freqMap] = op_CSIB0Correction(ftSpec_w);
% [ftSpec_B0corr] = op_CSIB0Correction(ftSpec_rmw,phaseMap,freqMap);
% 
% 
% %% === APODIZATION ===
% fprintf('Applying Gaussian apodization...\n');
% ftSpec_smooth = op_CSIApodize(ftSpec_B0corr, 'functionType', 'gaussian', 'fullWidthHalfMax', 25);
% ftSpec_smooth_w = op_CSIApodize(ftSpec_B0corr_w, 'functionType', 'gaussian', 'fullWidthHalfMax', 25);

%% === PLOT ===
op_CSIPlot(ftSpec_smooth);
% 
% fprintf('Generating brain and lipid masks...\n');
% [brain_area_raw, lipid_ring_raw] = createBrainArea(ftSpec_smooth, ...
%     ftSpec_smooth.sz(ftSpec_smooth.dims.x)/2, ...
%     ftSpec_smooth.sz(ftSpec_smooth.dims.y)/2);
% ftSpec_smooth.mask.brainmasks = brain_area_raw;
% ftSpec_smooth.mask.lipmasks = repmat(lipid_ring_raw, [1, 1, 5]);
% 
% io_CSIwritelcm(ftSpec_smooth_w,fullfile('C:\Users\divya\Downloads\fida2\new','lcm_ftSpec_smooth_w'));
% io_CSIwritelcm(ftSpec_smooth,fullfile('C:\Users\divya\Downloads\fida2\new','lcm_ftSpec_smooth'));

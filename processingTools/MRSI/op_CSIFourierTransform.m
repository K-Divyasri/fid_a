
function MRSIStruct = op_CSIFourierTransform(MRSIStruct, k_file, fourierTransform)
    arguments
        MRSIStruct (1,1) struct
        k_file (1,:) char {mustBeFileorDefault} = ""
        fourierTransform.spatial (1,1) logical {mustHaveSpatial(fourierTransform.spatial, MRSIStruct)}
        fourierTransform.spectral (1,1) logical {mustHaveSpectral(fourierTransform.spectral, MRSIStruct)}
    end

    % -- Match Brenden
    % 's defaulting logic exactly --
    fourierTransform = setDefaultFlags(fourierTransform, MRSIStruct);

    % --- Spatial Transform (if flagged) ---F
    if fourierTransform.spatial
        disp('Calculating spatial dimension');
        if (k_file == "")
            % Cartesian => fast Fourier Transform
            MRSIStruct = applyFastFourierTransformSpatial(MRSIStruct);
            Nt = MRSIStruct.sz(MRSIStruct.dims.t);
            MRSIStruct.spectralDwellTime = MRSIStruct.adcDwellTime;      % 1 time sample per ADC tick
            MRSIStruct.spectralWidth     = 1 / MRSIStruct.spectralDwellTime;
            MRSIStruct.spectralTime      = (0:Nt-1) * MRSIStruct.spectralDwellTime;
            MRSIStruct.adcTime      = (0:Nt-1) * MRSIStruct.adcDwellTime;
        else
            dComp=MRSIStruct;
            if ndims(dComp.data) == 5
                % --- 5D case ---
                % Original: [576 16 4 126 63]
                % Want: bring 126 to the front → permute to [126 576 16 4 63]
                % Then reshape 126*576 = 72576 → [72576 16 4 63]
                dComp.data = permute(dComp.data, [4 1 2 3 5]);
                dComp.data = reshape(dComp.data, [MRSIStruct.sz(4)*MRSIStruct.sz(1), MRSIStruct.sz(2), MRSIStruct.sz(3), MRSIStruct.sz(5)]);
                dComp.sz = size(dComp.data);
            
                % Update dims
                dComp.dims.t = 1;
                dComp.dims.coils = 2;
                dComp.dims.averages = 3;
                dComp.dims.ky = 4;
            
                % Zero out everything else
                fn = fieldnames(dComp.dims);
                for i = 1:numel(fn)
                    if ~ismember(fn{i}, {'t','coils','averages','ky'})
                        dComp.dims.(fn{i}) = 0;
                    end
                end
            
            elseif ndims(dComp.data) == 4
                % --- 4D case ---
                % Original: [576 16 126 63]
                % Want: bring 126 to the front → permute to [126 576 16 63]
                % Then reshape 126*576 = 72576 → [72576 16 63]
                dComp.data = permute(dComp.data, [3 1 2 4]);
                dComp.data = reshape(dComp.data, [MRSIStruct.sz(3)*MRSIStruct.sz(1), MRSIStruct.sz(2), MRSIStruct.sz(4)]);
                dComp.sz = size(dComp.data);
            
                % Update dims
                dComp.dims.t = 1;
                dComp.dims.coils = 2;
                dComp.dims.ky = 3;
            
                % Zero out everything else
                fn = fieldnames(dComp.dims);
                for i = 1:numel(fn)
                    if ~ismember(fn{i}, {'t','coils','ky'})
                        dComp.dims.(fn{i}) = 0;
                    end
                end
            
            else
                warning('Unexpected data dimensionality: %dD', ndims(dComp.data));
            end
            MRSIStruct=dComp;
            % Non-Cartesian => slow transform
            [kTable, kArray]  = readKFile(k_file);
            kPtsPerCycle      = getKPtsPerCycle(kTable);
            NPtemporal        = getTemporalPts(kTable, MRSIStruct);
            
            MRSIStruct = slowFourierTransform(MRSIStruct, kArray, kPtsPerCycle, NPtemporal);
            MRSIStruct = calculateSpectralValues(MRSIStruct, kPtsPerCycle, NPtemporal);
            Nt = MRSIStruct.sz(MRSIStruct.dims.t);
            MRSIStruct.adcTime      = (0:Nt-1) * MRSIStruct.adcDwellTime;
        end
        MRSIStruct = setFlags(MRSIStruct, 'spatialFT', true);
    end

    % --- Spectral Transform (if flagged) ---
    if fourierTransform.spectral
        disp('Calculating spectral dimension');
        MRSIStruct = fastFourierTransformTime_Fixed(MRSIStruct); %_fixed
    end
end

%% ===== MAIN FIX 1: Corrected calculateSpectralValues =====
function MRSIStruct = calculateSpectralValues_Fixed(MRSIStruct, kPtsPerCycle, NPtemporal)
    fprintf('=== calculateSpectralValues_Fixed ===\n');
    
    % CRITICAL FIX: For Rosette sequences, the spectral dwell time calculation
    % must account for the actual temporal sampling, not just k-space points
    
    adcDwellTime = MRSIStruct.adcDwellTime; % 5e-6 seconds
    
    % FIXED: For Rosette data, spectral dwell time = ADC dwell time * k-points per cycle
    % This gives the time between spectral samples after spatial reconstruction
    spectralDwellTime = adcDwellTime * kPtsPerCycle; % 5e-6 * 126 = 6.3e-4 s
    
    % Calculate spectral parameters
    spectralWidth = 1/spectralDwellTime; % Hz
    
    % CRITICAL FIX: Use the actual time dimension size (576), not derived values
    actualTimePoints = MRSIStruct.sz(MRSIStruct.dims.t); % Should be 576
    spectralTime = (0:(actualTimePoints-1)) * spectralDwellTime;
    
    fprintf('FIXED Parameters:\n');
    fprintf('  ADC dwell time: %.6e s\n', adcDwellTime);
    fprintf('  K-points per cycle: %d\n', kPtsPerCycle);
    fprintf('  Actual time points: %d\n', actualTimePoints);
    fprintf('  Spectral dwell time: %.6e s\n', spectralDwellTime);
    fprintf('  Spectral width: %.2f Hz\n', spectralWidth);
    fprintf('  Spectral time length: %d\n', length(spectralTime));
    
    % Set corrected values
    MRSIStruct.spectralWidth = spectralWidth;
    MRSIStruct.spectralDwellTime = spectralDwellTime;
    MRSIStruct.spectralTime = spectralTime;
    
    fprintf('=== calculateSpectralValues_Fixed END ===\n');
end

%% ===== MAIN FIX 2: Corrected op_NUFFTSpatial =====
function ftSpatial = op_NUFFTSpatial_Fixed(dComp, kFile_path)
    fprintf('=== op_NUFFTSpatial_Fixed START ===\n');
    
    % Load k-trajectory (keep existing logic)
    [kTable, ~] = readKFile_Fixed(kFile_path);
    
    % Handle different k-file formats
    if istable(kTable)
        if ismember('Kx', kTable.Properties.VariableNames) && ismember('Ky', kTable.Properties.VariableNames)
            kCoords = [kTable.Kx, kTable.Ky];
        else
            kCoords = [kTable{:, 2}, kTable{:, 3}];
        end
    else
        kCoords = kTable(:, 2:3);
    end
    
    % Normalize k-space coordinates
    kCoords = (kCoords / max(abs(kCoords(:)))) * pi;
    
    % Get dimensions
    sz = dComp.sz;
    dims = dComp.dims;
    nKshot = sz(dims.kshot);
    nKpts = sz(dims.kpts);
    totalKPoints = nKshot * nKpts;
    
    fprintf('Data dimensions: %s\n', mat2str(sz));
    fprintf('K-space: %d shots × %d points = %d total\n', nKshot, nKpts, totalKPoints);
    
    % Verify k-space dimensions
    if size(kCoords, 1) ~= totalKPoints
        error('K-space file length (%d) ≠ data k-space dimensions (%d)', ...
              size(kCoords, 1), totalKPoints);
    end
    
    % NUFFT setup
    Nx = length(dComp.coordinates.x);
    Ny = length(dComp.coordinates.y);
    
    st = nufft_init(kCoords, [Nx, Ny], [6, 6], [2*Nx, 2*Ny], [Nx, Ny]/2);
    
    % Process data based on dimensions
    hasAverages = (dims.averages > 0);
    
    if hasAverages && ndims(dComp.data) == 5
        % 5D: [t, coils, averages, kpts, kshot]
        data_reshaped = reshape(dComp.data, [sz(1)*sz(2)*sz(3), sz(4)*sz(5)]);
        data_transposed = data_reshaped.';
        
        imgFlat = nufft_adj(data_transposed, st);
        imgData = reshape(imgFlat, [Nx, Ny, sz(1), sz(2), sz(3)]);
        imgData = permute(imgData, [3, 4, 5, 2, 1]); % [t, coils, averages, x, y]
        
        new_dims = dims;
        new_dims.x = 5;
        new_dims.y = 4;
        new_dims.kshot = 0;
        new_dims.kpts = 0;
        
    else % 4D: [t, coils, kpts, kshot]
        data_reshaped = reshape(dComp.data, [sz(1)*sz(2), sz(3)*sz(4)]);
        data_transposed = data_reshaped.';
        
        imgFlat = nufft_adj(data_transposed, st);
        imgData = reshape(imgFlat, [Nx, Ny, sz(1), sz(2)]);
        imgData = permute(imgData, [3, 4, 2, 1]); % [t, coils, x, y]
        
        new_dims = dims;
        new_dims.x = 4;
        new_dims.y = 3;
        new_dims.kshot = 0;
        new_dims.kpts = 0;
    end
    
    % Update structure
    ftSpatial = dComp;
    ftSpatial.data = imgData;
    ftSpatial.sz = size(imgData);
    ftSpatial.dims = new_dims;
    ftSpatial.flags.spatialFT = 1;
    
    % FIXED: Spectral parameter calculation
    kPtsPerCycle = getKPtsPerCycle_Fixed(kTable);
    NPtemporal = sz(dims.t); % Use actual time dimension size
    
    ftSpatial = calculateSpectralValues_Fixed(ftSpatial, kPtsPerCycle, NPtemporal);
    
    fprintf('=== op_NUFFTSpatial_Fixed END ===\n');
end

%% ===== MAIN FIX 3: Corrected Spectral Fourier Transform =====
function MRSIStruct = op_CSIFourierTransform_Fixed(MRSIStruct, varargin)
    % Fixed spectral Fourier transform with proper apodization
    
    p = inputParser;
    addRequired(p, 'MRSIStruct');
    addParameter(p, 'apodization', 3, @isnumeric); % Line broadening in Hz
    addParameter(p, 'zeroFill', 1, @isnumeric);    % Zero-filling factor
    parse(p, MRSIStruct, varargin{:});
    
    fprintf('=== op_CSIFourierTransform_Fixed ===\n');
    
    data = MRSIStruct.data;
    timeDim = getDimension(MRSIStruct, 't');
    
    if timeDim == 0
        error('No time dimension found for spectral FT');
    end
    
    fprintf('Input data size: %s\n', mat2str(size(data)));
    fprintf('Time dimension: %d\n', timeDim);
    
    % STEP 1: Apply Gaussian apodization to reduce ripples
    if p.Results.apodization > 0
        data = applyGaussianApodization_Fixed(data, timeDim, MRSIStruct.adcDwellTime, p.Results.apodization);
        fprintf('Applied %.1f Hz Gaussian apodization\n', p.Results.apodization);
    end
    
    % STEP 2: Optional zero-filling for better spectral resolution
    if p.Results.zeroFill > 1
        data = applyZeroFilling_Fixed(data, timeDim, p.Results.zeroFill);
        fprintf('Applied %dx zero-filling\n', p.Results.zeroFill);
    end
    
    % STEP 3: Spectral Fourier Transform
    % For MRS: IFFT transforms time domain to frequency domain
    data_ft = ifft(data, [], timeDim);
    
    % STEP 4: FFT shift to center the spectrum
    data_ft = fftshift(data_ft, timeDim);
    
    % STEP 5: Apply automatic phase correction
    data_ft = applyAutoPhaseCorrection_Fixed(data_ft, timeDim);
    
    % Update structure
    MRSIStruct.data = data_ft;
    MRSIStruct.sz = size(data_ft);
    
    % STEP 6: Calculate corrected PPM axis
    ppm = calculatePPM_Fixed(MRSIStruct);
    MRSIStruct.ppm = ppm;
    
    % Update flags and dimensions
    MRSIStruct.flags.spectralFT = 1;
    
    % Rename dimension from 't' to 'f'
    fDim = getDimension(MRSIStruct, 't');
    MRSIStruct.dims.f = fDim;
    MRSIStruct.dims.t = 0;
    
    fprintf('Spectral FT completed. Output size: %s\n', mat2str(size(data_ft)));
    fprintf('PPM range: %.2f to %.2f ppm\n', min(ppm), max(ppm));
    fprintf('=== op_CSIFourierTransform_Fixed END ===\n');
end

%% ===== SUPPORTING FUNCTIONS =====

function data_apod = applyGaussianApodization_Fixed(data, timeDim, dwellTime, lbHz)
    % Apply Gaussian line broadening to reduce Gibbs ringing
    
    nPts = size(data, timeDim);
    t = (0:(nPts-1)) * dwellTime;
    
    % Gaussian apodization function: exp(-π * lb * t)
    gaussFunc = exp(-pi * lbHz * t);
    
    % Reshape for broadcasting
    shape = ones(1, ndims(data));
    shape(timeDim) = length(gaussFunc);
    gaussFunc = reshape(gaussFunc, shape);
    
    data_apod = data .* gaussFunc;
end

function data_zf = applyZeroFilling_Fixed(data, timeDim, zfFactor)
    % Apply zero-filling in time domain
    
    sz = size(data);
    sz_new = sz;
    sz_new(timeDim) = sz_new(timeDim) * zfFactor;
    
    data_zf = zeros(sz_new, 'like', data);
    
    % Copy original data to beginning of zero-filled array
    idx = cell(1, ndims(data));
    for i = 1:ndims(data)
        if i == timeDim
            idx{i} = 1:sz(i);
        else
            idx{i} = 1:sz(i);
        end
    end
    
    data_zf(idx{:}) = data;
end

function data_corrected = applyAutoPhaseCorrection_Fixed(data, timeDim)
    % Apply automatic zero-order phase correction
    
    % Find peak in average spectrum
    avgSpec = squeeze(mean(data, [2:ndims(data)]));
    if timeDim ~= 1
        avgSpec = squeeze(mean(data, setdiff(1:ndims(data), timeDim)));
    end
    
    [~, peakIdx] = max(abs(avgSpec));
    
    % Get phase at peak
    if timeDim == 1
        peakPhase = angle(data(peakIdx, :));
    else
        peakData = data;
        for dim = 1:ndims(data)
            if dim ~= timeDim
                peakData = mean(peakData, dim);
            end
        end
        peakPhase = angle(peakData);
    end
    
    % Calculate average phase
    meanPhase = angle(mean(exp(1i * peakPhase(:))));
    
    % Apply phase correction
    data_corrected = data * exp(-1i * meanPhase);
    
    fprintf('Applied phase correction: %.1f degrees\n', meanPhase * 180/pi);
end

function ppm = calculatePPM_Fixed(MRSIStruct)
    % Calculate corrected PPM axis
    
    spectralWidth = MRSIStruct.spectralWidth; % Hz
    nPoints = MRSIStruct.sz(MRSIStruct.dims.f);
    Bo = MRSIStruct.Bo; % Tesla
    gamma = MRSIStruct.gamma; % MHz/T
    
    % Create frequency axis (centered at 0)
    df = spectralWidth / nPoints;
    freqs = (-spectralWidth/2 : df : spectralWidth/2 - df);
    
    % Convert to PPM
    ppm = freqs / (gamma * Bo * 1e6);
    
    % Add water reference (4.7 ppm for 1H)
    if strcmp(MRSIStruct.nucleus, '1H')
        ppm = ppm + 4.7;
    end
end

function kPtsPerCycle = getKPtsPerCycle_Fixed(kTable)
    % Extract k-points per cycle with better error handling
    
    try
        if istable(kTable)
            if ismember('TR', kTable.Properties.VariableNames)
                num_TR = max(kTable.TR);
                kPtsPerCycle = height(kTable) / num_TR;
            else
                kPtsPerCycle = 126; % Default for Rosette
            end
        elseif ismatrix(kTable) && size(kTable, 2) >= 5
            num_TR = max(kTable(:, 5));
            if num_TR > 0
                kPtsPerCycle = size(kTable, 1) / num_TR;
            else
                kPtsPerCycle = 126;
            end
        else
            kPtsPerCycle = 126; % Default
        end
    catch
        kPtsPerCycle = 126; % Fallback
    end
    
    fprintf('K-points per cycle: %d\n', kPtsPerCycle);
end

function [kTable, kArray] = readKFile_Fixed(kFileName)
    % Read k-space file with robust error handling
    
    if isempty(kFileName) || ~isfile(kFileName)
        kTable = [];
        kArray = [];
        return;
    end
    
    try
        kTable = readtable(kFileName);
        if ismember('Kx', kTable.Properties.VariableNames) && ismember('Ky', kTable.Properties.VariableNames)
            kArray = [kTable.Kx, kTable.Ky];
        else
            kArray = [kTable{:, 2}, kTable{:, 3}];
        end
    catch
        try
            kTable = readmatrix(kFileName);
            kArray = kTable(:, 2:3);
        catch ME
            warning('Failed to read k-file: %s', ME.message);
            kTable = [];
            kArray = [];
        end
    end
end

%% ===== HELPER FUNCTIONS (from your existing code) =====

function dim = getDimension(MRSIStruct, dimName)
    if isfield(MRSIStruct.dims, dimName)
        dim = MRSIStruct.dims.(dimName);
    else
        dim = 0;
    end
end

function adcDwellTime = getAdcDwellTime(MRSIStruct)
    adcDwellTime = MRSIStruct.adcDwellTime;
end

%% ===== COMPLETE CORRECTED PROCESSING PIPELINE =====
function process_corrected_pipeline()
    % Example of how to use the corrected functions
    
    fprintf('=== CORRECTED PROCESSING PIPELINE ===\n');
    
    % Your existing variables (from the debugging output):
    % dComp: [576×16×6×126×76] - metabolite data with averages
    % dComp_w: [576×16×126×76] - water reference without averages
    % kFile: path to k-space trajectory file
    
    % STEP 1: Apply corrected NUFFT spatial transform
    fprintf('Step 1: Spatial NUFFT...\n');
    ftSpatial = op_NUFFTSpatial_Fixed(dComp, kFile);
    ftSpatial_w = op_NUFFTSpatial_Fixed(dComp_w, kFile);
    
    % STEP 2: Coil combination (keep your existing function)
    fprintf('Step 2: Coil combination...\n');
    [coilCombined_w, phase, weights] = op_CSICombineCoils1(ftSpatial_w);
    coilCombined = op_CSICombineCoils1(ftSpatial, 1, phase, weights);
    
    % STEP 3: Average combination (keep your existing function)
    fprintf('Step 3: Averaging...\n');
    ccav = op_CSIAverage(coilCombined);
    ccav_w = op_CSIAverage(coilCombined_w);
    
    % STEP 4: Apply corrected spectral Fourier transform
    fprintf('Step 4: Spectral FT with corrections...\n');
    ftSpec = op_CSIFourierTransform_Fixed(ccav, 'apodization', 3, 'zeroFill', 1);
    ftSpec_w = op_CSIFourierTransform_Fixed(ccav_w, 'apodization', 3, 'zeroFill', 1);
    
    % STEP 5: Additional processing (lipid removal, B0 correction, etc.)
    fprintf('Step 5: Additional processing...\n');
    % Continue with your existing functions...
    
    % STEP 6: Final apodization for display
    fprintf('Step 6: Final smoothing...\n');
    ftSpec_smooth = op_CSIApodize(ftSpec, 'functionType', 'gaussian', 'fullWidthHalfMax', 15);
    
    % STEP 7: Plot results
    fprintf('Step 7: Plotting...\n');
    op_CSIPlot(ftSpec_smooth);
    
    fprintf('=== CORRECTED PIPELINE COMPLETE ===\n');
end

%% ===== QUICK FIX FOR YOUR CURRENT PIPELINE =====
function quick_fix_current_data()
    % If you want to quickly fix your current ftSpec data
    
    % Load your current ftSpec variable
    % ftSpec = ... (your current data)
    
    fprintf('=== QUICK FIX FOR CURRENT DATA ===\n');
    
    % Recalculate spectral parameters
    ftSpec.spectralDwellTime = ftSpec.adcDwellTime * 126; % 5e-6 * 126
    ftSpec.spectralWidth = 1/ftSpec.spectralDwellTime;
    
    % Recalculate spectral time
    nPts = ftSpec.sz(1); % 576
    ftSpec.spectralTime = (0:(nPts-1)) * ftSpec.spectralDwellTime;
    
    % Recalculate PPM axis
    df = ftSpec.spectralWidth / nPts;
    freqs = (-ftSpec.spectralWidth/2 : df : ftSpec.spectralWidth/2 - df);
    ftSpec.ppm = freqs / (ftSpec.gamma * ftSpec.Bo * 1e6) + 4.7;
    
    % Apply additional smoothing to reduce remaining ripples
    ftSpec_smooth = op_CSIApodize(ftSpec, 'functionType', 'gaussian', 'fullWidthHalfMax', 10);
    
    % Plot corrected result
    op_CSIPlot(ftSpec_smooth);
    
    fprintf('Quick fix applied. Try plotting ftSpec_smooth.\n');
    fprintf('=== QUICK FIX COMPLETE ===\n');
end

%% ===== DIAGNOSTIC FUNCTION =====
function diagnose_spectral_parameters(MRSIStruct)
    % Diagnostic function to check spectral parameters
    
    fprintf('=== SPECTRAL PARAMETERS DIAGNOSIS ===\n');
    fprintf('Data size: %s\n', mat2str(MRSIStruct.sz));
    fprintf('ADC dwell time: %.6e s\n', MRSIStruct.adcDwellTime);
    fprintf('Spectral dwell time: %.6e s\n', MRSIStruct.spectralDwellTime);
    fprintf('Spectral width: %.2f Hz\n', MRSIStruct.spectralWidth);
    fprintf('Spectral time length: %d\n', length(MRSIStruct.spectralTime));
    fprintf('PPM range: %.2f to %.2f ppm\n', min(MRSIStruct.ppm), max(MRSIStruct.ppm));
    
    % Check for common issues
    expectedSpectralDT = MRSIStruct.adcDwellTime * 126;
    if abs(MRSIStruct.spectralDwellTime - expectedSpectralDT) > 1e-8
        fprintf('WARNING: Spectral dwell time seems incorrect!\n');
        fprintf('Expected: %.6e s, Got: %.6e s\n', expectedSpectralDT, MRSIStruct.spectralDwellTime);
    end
    
    fprintf('=== DIAGNOSIS COMPLETE ===\n');
end




%% -------------------------------------------------------------------------
% sft2_Operator: kinda, but vectorized
%% -------------------------------------------------------------------------
function sft2_Oper = sft2_Operator(InTraj, OutTraj, Ift_flag)
    % In the original code:
    %   if (~Ift_flag) => exponent = -2*pi*1i
    %   else exponent =  2*pi*1i
    %   if Ift_flag => divide by Nx as well

    if ~Ift_flag
        Expy = -2*pi*1i;    % forward transform exponent
    else
        Expy =  2*pi*1i;    % inverse transform exponent
    end

    NOut   = size(OutTraj,1);
    NIn    = size(InTraj,1);

    % Vectorized approach
    xTerm   = OutTraj(:,1)*InTraj(:,1)';  % size [NOut x NIn]
    yTerm   = OutTraj(:,2)*InTraj(:,2)';  % size [NOut x NIn]
    sft2_Oper = exp(Expy*(xTerm + yTerm));

    if Ift_flag
        sft2_Oper = sft2_Oper / NIn;  % divide by number of input k-points
    end
end


%% -------------------------------------------------------------------------
% mustHaveSpatial / mustHaveSpectral
%% -------------------------------------------------------------------------
function mustHaveSpatial(a, in)
    if isfield(in, 'spatialFT') && (a == true && in.spatialFT == 1)
        error('Spatial Fourier Transform already done!');
    end
end

function mustHaveSpectral(a, in)
    if isfield(in, 'spectral') && (a == true && in.spectral == 1)
        error('Spectral Fourier Transform already done!');
    end
end

function mustBeFileorDefault(file)
    if ~isfile(file) && ~strcmp(file, "")
        error('Invalid k_file, must be an existing file or empty.');
    end
end


%% -------------------------------------------------------------------------
% setDefaultFlags: same stuff
%% -------------------------------------------------------------------------
function fourierTransform = setDefaultFlags(fourierTransform, in)
    if ~isfield(fourierTransform, 'spatial')
        if in.flags.spatialFT
            fourierTransform.spatial = 0;
        else
            fourierTransform.spatial = 1;
        end
    end
    if ~isfield(fourierTransform, 'spectral')
        if in.flags.spectralFT
            fourierTransform.spectral = 0;
        else
            fourierTransform.spectral = 1;
        end
    end
end


%% -------------------------------------------------------------------------
% applyFastFourierTransformSpatial: same stuff
%% -------------------------------------------------------------------------
function MRSIStruct = applyFastFourierTransformSpatial(MRSIStruct)
    disp('Applying fast fourier transform (Cartesian)');

    % 1) half-pixel shift
    MRSIStruct = halfPixelShift(MRSIStruct);

    % 2) FFT along x
    data = getData(MRSIStruct);
    xDim = getDimension(MRSIStruct, 'kx');
    if mod(getSizeFromDimensions(MRSIStruct, {'kx'}), 2) == 1
        data = circshift(data, 1, xDim);
    end
    data = fftshift( fft( fftshift(data, xDim), [], xDim ), xDim);

    % 3) FFT along y
    yDim = getDimension(MRSIStruct, 'ky');
    if mod(getSizeFromDimensions(MRSIStruct, {'ky'}), 2) == 1
        data = circshift(data, 1, yDim);
    end
    data = fftshift( fft( fftshift(data, yDim), [], yDim ), yDim);

    MRSIStruct = setData(MRSIStruct, data);

    % 4) re-label dims
    MRSIStruct = setDimension(MRSIStruct, 'x',  getDimension(MRSIStruct, 'kx'));
    MRSIStruct = setDimension(MRSIStruct, 'y',  getDimension(MRSIStruct, 'ky'));
    MRSIStruct = setDimension(MRSIStruct, 'z',  getDimension(MRSIStruct, 'kz'));
    MRSIStruct = setDimension(MRSIStruct, 'kx', 0);
    MRSIStruct = setDimension(MRSIStruct, 'ky', 0);
    MRSIStruct = setDimension(MRSIStruct, 'kz', 0);
end


%% -------------------------------------------------------------------------
% slowFourierTransform: same chunk approach as Brenden, preserving dims
%% -------------------------------------------------------------------------
function MRSIStruct = slowFourierTransform(MRSIStruct, kTrajectory, kPtsPerCycle, NPtemporal)
    [xCoordinates, yCoordinates, imageTrajectory] = getImageTrajectory(MRSIStruct);

    % Build slow transform operator
    sftOperator = sft2_Operator(kTrajectory, imageTrajectory, 1);
    [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'t','ky'});

    data = getData(MRSIStruct);

    % Apply chunked slow transform
    image = applySlowFourierTranformMatrix(MRSIStruct, sftOperator, data, NPtemporal, kPtsPerCycle);
    MRSIStruct = setData(MRSIStruct, image);

    % Re-label dims
    kyDimension = getDimension(MRSIStruct, 'ky');
    prevPermute = removeDimPrevPermute(prevPermute, kyDimension);
    prevPermute = addDimPrevPermute(prevPermute, 'y', kyDimension);
    prevPermute = addDimPrevPermute(prevPermute, 'x', kyDimension + 1);

    % Match Brenden's logic for final size
    prevSize(1) = NPtemporal;
    prevSize(2) = length(yCoordinates);
    prevSize    = [prevSize(1:2), length(xCoordinates), prevSize(3:end)];

    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);
end


%% -------------------------------------------------------------------------
% applySlowFourierTranformMatrix: chunk over time blocks, same stuff
%% -------------------------------------------------------------------------
function outData = applySlowFourierTranformMatrix(MRSIStruct, sftOperator, data, NPtemporal, kPtsPerCycle)
    yLength = length(getCoordinates(MRSIStruct, 'y'));
    xLength = length(getCoordinates(MRSIStruct, 'x'));

    extrasSize = getSizeFromDimensions(MRSIStruct, {'extras'});
    outDims    = [NPtemporal, yLength, xLength, extrasSize];

    outData = zeros(outDims);  % same numeric class as default double

    for iPoint = 1:NPtemporal
        startPt = (iPoint - 1)*kPtsPerCycle + 1;
        endPt   = iPoint * kPtsPerCycle;

        kSlice = data(startPt:endPt, :, :);
        vectSlice = reshape(kSlice, [], size(kSlice,3));  % flatten

        ftSlice = sftOperator * vectSlice;
        ftSlice = reshape(ftSlice, [yLength, xLength, size(ftSlice,2)]);

        outData(iPoint, :, :, :) = ftSlice;
    end
end


%% -------------------------------------------------------------------------
% getImageTrajectory
%% -------------------------------------------------------------------------
function [xCoordinates, yCoordinates, imageTrajectory] = getImageTrajectory(MRSIStruct)
    xCoordinates = getCoordinates(MRSIStruct, 'x');
    yCoordinates = getCoordinates(MRSIStruct, 'y');

    [xx, yy]      = meshgrid(xCoordinates, yCoordinates);
    imageTrajectory = [xx(:), yy(:)];
end


%% -------------------------------------------------------------------------
% fastFourierTransformTime: same stuff
%% -------------------------------------------------------------------------
function MRSIStruct = fastFourierTransformTime_Fixed(MRSIStruct)
    data = getData(MRSIStruct);
    timeDimension = getDimension(MRSIStruct, 't');
    
    % Fourier transform in the spectral domain (time -> frequency)
    data = fftshift(ifft(data, [], timeDimension), timeDimension);
    MRSIStruct = setData(MRSIStruct, data);

    % Compute PPM axis
    ppm = calculatePPM(MRSIStruct);
    if strcmp(MRSIStruct.nucleus,'1H')
        ppm = ppm + 4.65;
    end
    MRSIStruct = setPPM(MRSIStruct, ppm);

    % Update flags
    MRSIStruct = setFlags(MRSIStruct, 'spectralFT', true);

    % ---- NEW: Rename dimension label from 't' to 'f' ----
    fDim = getDimension(MRSIStruct, 't');  % Get dimension index for 't'
    MRSIStruct = setDimension(MRSIStruct, 'f', fDim);  % Set same index for 'f'
    MRSIStruct = setDimension(MRSIStruct, 't', 0);      % Clear old 't' label
end


%% -------------------------------------------------------------------------
% calculatePPM: same stuff
%% -------------------------------------------------------------------------
function ppmVals = calculatePPM(MRSIStruct)
    gammaVal      = MRSIStruct.gamma;
    spectralWidth = getSpectralWidth(MRSIStruct);
    timeSize      = getSizeFromDimensions(MRSIStruct, {'t'});

    step       = spectralWidth / timeSize;
    lowerBound = -spectralWidth/2 + step/2;
    upperBound =  spectralWidth/2 - step/2;

    frequencyArray = lowerBound : step : upperBound;
    ppmVals        = -frequencyArray / (MRSIStruct.Bo * gammaVal);
end


%% -------------------------------------------------------------------------
% calculateSpectralValues: same stuff
%% -------------------------------------------------------------------------
% function MRSIStruct = calculateSpectralValues(MRSIStruct, kPtsPerCycle, NPtemporal)
%     spectralDwellTime = calculateSpectralDwellTime(MRSIStruct, kPtsPerCycle);
%     spectralWidth     = 1 / spectralDwellTime;
%     spectralTime      = calculateSpectralTime(spectralDwellTime, NPtemporal);
% 
%     MRSIStruct = setSpectralWidth(MRSIStruct, spectralWidth);
%     MRSIStruct = setSpectralDwellTime(MRSIStruct, spectralDwellTime);
%     MRSIStruct = setSpectralTime(MRSIStruct, spectralTime);
% end
function MRSIStruct = calculateSpectralValues(MRSIStruct, kPtsPerCycle, NPtemporal)
    fprintf('=== calculateSpectralValues (FIXED) ===\n');
    
    % CRITICAL FIX: Use NPtemporal for spectral calculations, not kPtsPerCycle
    % The spectral dwell time should be based on temporal points, not k-space points
    adcDwellTime = MRSIStruct.adcDwellTime;
    
    % FIXED: Spectral dwell time calculation
    % For Rosette data: spectral_dt = adc_dt * k_points_per_temporal_cycle
    spectralDwellTime = adcDwellTime * kPtsPerCycle;
    
    % Calculate spectral width and time vector
    spectralWidth = 1/spectralDwellTime;
    spectralTime = (0:(NPtemporal-1)) * spectralDwellTime;
    
    fprintf('Debug (FIXED): adcDwellTime = %.6e s\n', adcDwellTime);
    fprintf('Debug (FIXED): kPtsPerCycle = %d\n', kPtsPerCycle);
    fprintf('Debug (FIXED): NPtemporal = %d\n', NPtemporal);
    fprintf('Debug (FIXED): spectralDwellTime = %.6e s\n', spectralDwellTime);
    fprintf('Debug (FIXED): spectralWidth = %.2f Hz\n', spectralWidth);
    fprintf('Debug (FIXED): spectralTime length = %d\n', length(spectralTime));
    
    % Set values in structure
    MRSIStruct.spectralWidth = spectralWidth;
    MRSIStruct.spectralDwellTime = spectralDwellTime;
    MRSIStruct.spectralTime = spectralTime;
    
    fprintf('=== calculateSpectralValues (FIXED) END ===\n');
end


function spectralDwellTime = calculateSpectralDwellTime(MRSIStruct, spatialPoints)
    adcDwellTime = getAdcDwellTime(MRSIStruct);
    spectralDwellTime = spatialPoints * adcDwellTime;
end

function tAxis = calculateSpectralTime(spectralDwellTime, spatialPoints)
    tAxis = 0:spectralDwellTime:spectralDwellTime*(spatialPoints - 1);
end


%% -------------------------------------------------------------------------
% halfPixelShift: same stuff
%% -------------------------------------------------------------------------
function MRSIStruct = halfPixelShift(MRSIStruct)
    kx = getCoordinates(MRSIStruct, 'kx');
    ky = getCoordinates(MRSIStruct, 'ky');

    halfPixelX = getVoxSize(MRSIStruct, 'x')/2;
    halfPixelY = getVoxSize(MRSIStruct, 'y')/2;

    kShift = kx*halfPixelX + ky'*halfPixelY;

    [MRSIStruct, prevPerm, prevSz] = reshapeDimensions(MRSIStruct, {'ky','kx'});
    data = getData(MRSIStruct);

    data = data .* exp(-1i*2*pi*kShift);

    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prevPerm, prevSz);
end







%% ===== FIX 6: Additional helper functions =====

% Zero-filling function (optional enhancement)
function data_zf = zero_fill_spectral(data, timeDim, zf_factor)
    % Zero-fill in spectral dimension to improve resolution
    sz = size(data);
    sz_new = sz;
    sz_new(timeDim) = sz_new(timeDim) * zf_factor;
    
    data_zf = zeros(sz_new, 'like', data);
    
    % Copy original data to center of zero-filled array
    idx = cell(1, ndims(data));
    for i = 1:ndims(data)
        if i == timeDim
            start_idx = floor((sz_new(i) - sz(i))/2) + 1;
            idx{i} = start_idx:(start_idx + sz(i) - 1);
        else
            idx{i} = 1:sz(i);
        end
    end
    
    data_zf(idx{:}) = data;
end
% % op_CSIFourierTransform.m
% % Brenden Kadota, SunnyBrook Hospital, 2021 - vectorised
% %
% % Performs spatial & spectral Fourier transforms on MRSI data.
% %   - If k_file is empty => Cartesian => fast FFT on (kx, ky).
% %   - If k_file is provided => non-Cartesian => chunked slow transform.
% %   - Spectral transform => IFFT along 'time' dimension.
% 
% function MRSIStruct = op_CSIFourierTransform(MRSIStruct, k_file, fourierTransform)
%     arguments
%         MRSIStruct (1,1) struct
%         k_file (1,:) char {mustBeFileorDefault} = ""
%         fourierTransform.spatial (1,1) logical {mustHaveSpatial(fourierTransform.spatial, MRSIStruct)}
%         fourierTransform.spectral (1,1) logical {mustHaveSpectral(fourierTransform.spectral, MRSIStruct)}
%     end
% 
%     % -- Match Brenden's defaulting logic exactly --
%     fourierTransform = setDefaultFlags(fourierTransform, MRSIStruct);
% 
%     % --- Spatial Transform (if flagged) ---
%     if fourierTransform.spatial
%         disp('Calculating spatial dimension');
%         if (k_file == "")
%             % Cartesian => fast Fourier Transform
%             MRSIStruct = applyFastFourierTransformSpatial(MRSIStruct);
%         else
%             % Non-Cartesian => slow transform
%             [kTable, kArray]  = readKFile(k_file);
%             kPtsPerCycle      = getKPtsPerCycle(kTable);
%             NPtemporal        = getTemporalPts(kTable, MRSIStruct);
% 
%             MRSIStruct = slowFourierTransform(MRSIStruct, kArray, kPtsPerCycle, NPtemporal);
%             MRSIStruct = calculateSpectralValues(MRSIStruct, kPtsPerCycle, NPtemporal);
%         end
%         MRSIStruct = setFlags(MRSIStruct, 'spatialFT', true);
%     end
% 
%     % --- Spectral Transform (if flagged) ---
%     if fourierTransform.spectral
%         disp('Calculating spectral dimension');
%         MRSIStruct = fastFourierTransformTime(MRSIStruct);
%     end
% end
% 
% 
% %% -------------------------------------------------------------------------
% % sft2_Operator: kinda, but vectorized
% %% -------------------------------------------------------------------------
% function sft2_Oper = sft2_Operator(InTraj, OutTraj, Ift_flag)
%     % In the original code:
%     %   if (~Ift_flag) => exponent = -2*pi*1i
%     %   else exponent =  2*pi*1i
%     %   if Ift_flag => divide by Nx as well
% 
%     if ~Ift_flag
%         Expy = -2*pi*1i;    % forward transform exponent
%     else
%         Expy =  2*pi*1i;    % inverse transform exponent
%     end
% 
%     NOut   = size(OutTraj,1);
%     NIn    = size(InTraj,1);
% 
%     % Vectorized approach
%     xTerm   = OutTraj(:,1)*InTraj(:,1)';  % size [NOut x NIn]
%     yTerm   = OutTraj(:,2)*InTraj(:,2)';  % size [NOut x NIn]
%     sft2_Oper = exp(Expy*(xTerm + yTerm));
% 
%     if Ift_flag
%         sft2_Oper = sft2_Oper / NIn;  % divide by number of input k-points
%     end
% end
% 
% 
% %% -------------------------------------------------------------------------
% % mustHaveSpatial / mustHaveSpectral
% %% -------------------------------------------------------------------------
% function mustHaveSpatial(a, in)
%     if isfield(in, 'spatialFT') && (a == true && in.spatialFT == 1)
%         error('Spatial Fourier Transform already done!');
%     end
% end
% 
% function mustHaveSpectral(a, in)
%     if isfield(in, 'spectral') && (a == true && in.spectral == 1)
%         error('Spectral Fourier Transform already done!');
%     end
% end
% 
% function mustBeFileorDefault(file)
%     if ~isfile(file) && ~strcmp(file, "")
%         error('Invalid k_file, must be an existing file or empty.');
%     end
% end
% 
% 
% %% -------------------------------------------------------------------------
% % setDefaultFlags: same stuff
% %% -------------------------------------------------------------------------
% function fourierTransform = setDefaultFlags(fourierTransform, in)
%     if ~isfield(fourierTransform, 'spatial')
%         if in.flags.spatialFT
%             fourierTransform.spatial = 0;
%         else
%             fourierTransform.spatial = 1;
%         end
%     end
%     if ~isfield(fourierTransform, 'spectral')
%         if in.flags.spectralFT
%             fourierTransform.spectral = 0;
%         else
%             fourierTransform.spectral = 1;
%         end
%     end
% end
% 
% 
% %% -------------------------------------------------------------------------
% % applyFastFourierTransformSpatial: same stuff
% %% -------------------------------------------------------------------------
% function MRSIStruct = applyFastFourierTransformSpatial(MRSIStruct)
%     disp('Applying fast fourier transform (Cartesian)');
% 
%     % 1) half-pixel shift
%     MRSIStruct = halfPixelShift(MRSIStruct);
% 
%     % 2) FFT along x
%     data = getData(MRSIStruct);
%     xDim = getDimension(MRSIStruct, 'kx');
%     if mod(getSizeFromDimensions(MRSIStruct, {'kx'}), 2) == 1
%         data = circshift(data, 1, xDim);
%     end
%     data = fftshift( fft( fftshift(data, xDim), [], xDim ), xDim);
% 
%     % 3) FFT along y
%     yDim = getDimension(MRSIStruct, 'ky');
%     if mod(getSizeFromDimensions(MRSIStruct, {'ky'}), 2) == 1
%         data = circshift(data, 1, yDim);
%     end
%     data = fftshift( fft( fftshift(data, yDim), [], yDim ), yDim);
% 
%     MRSIStruct = setData(MRSIStruct, data);
% 
%     % 4) re-label dims
%     MRSIStruct = setDimension(MRSIStruct, 'x',  getDimension(MRSIStruct, 'kx'));
%     MRSIStruct = setDimension(MRSIStruct, 'y',  getDimension(MRSIStruct, 'ky'));
%     MRSIStruct = setDimension(MRSIStruct, 'z',  getDimension(MRSIStruct, 'kz'));
%     MRSIStruct = setDimension(MRSIStruct, 'kx', 0);
%     MRSIStruct = setDimension(MRSIStruct, 'ky', 0);
%     MRSIStruct = setDimension(MRSIStruct, 'kz', 0);
% end
% 
% 
% %% -------------------------------------------------------------------------
% % slowFourierTransform: same chunk approach as Brenden, preserving dims
% %% -------------------------------------------------------------------------
% function MRSIStruct = slowFourierTransform(MRSIStruct, kTrajectory, kPtsPerCycle, NPtemporal)
%     [xCoordinates, yCoordinates, imageTrajectory] = getImageTrajectory(MRSIStruct);
% 
%     % Build slow transform operator
%     sftOperator = sft2_Operator(kTrajectory, imageTrajectory, 1);
%     [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'t','ky'});
% 
%     data = getData(MRSIStruct);
% 
%     % Apply chunked slow transform
%     image = applySlowFourierTranformMatrix(MRSIStruct, sftOperator, data, NPtemporal, kPtsPerCycle);
%     MRSIStruct = setData(MRSIStruct, image);
% 
%     % Re-label dims
%     kyDimension = getDimension(MRSIStruct, 'ky');
%     prevPermute = removeDimPrevPermute(prevPermute, kyDimension);
%     prevPermute = addDimPrevPermute(prevPermute, 'y', kyDimension);
%     prevPermute = addDimPrevPermute(prevPermute, 'x', kyDimension + 1);
% 
%     % Match Brenden's logic for final size
%     prevSize(1) = NPtemporal;
%     prevSize(2) = length(yCoordinates);
%     prevSize    = [prevSize(1:2), length(xCoordinates), prevSize(3:end)];
% 
%     MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);
% end
% 
% 
% %% -------------------------------------------------------------------------
% % applySlowFourierTranformMatrix: chunk over time blocks, same stuff
% %% -------------------------------------------------------------------------
% function outData = applySlowFourierTranformMatrix(MRSIStruct, sftOperator, data, NPtemporal, kPtsPerCycle)
%     yLength = length(getCoordinates(MRSIStruct, 'y'));
%     xLength = length(getCoordinates(MRSIStruct, 'x'));
% 
%     extrasSize = getSizeFromDimensions(MRSIStruct, {'extras'});
%     outDims    = [NPtemporal, yLength, xLength, extrasSize];
% 
%     outData = zeros(outDims);  % same numeric class as default double
% 
%     for iPoint = 1:NPtemporal
%         startPt = (iPoint - 1)*kPtsPerCycle + 1;
%         endPt   = iPoint * kPtsPerCycle;
% 
%         kSlice = data(startPt:endPt, :, :);
%         vectSlice = reshape(kSlice, [], size(kSlice,3));  % flatten
% 
%         ftSlice = sftOperator * vectSlice;
%         ftSlice = reshape(ftSlice, [yLength, xLength, size(ftSlice,2)]);
% 
%         outData(iPoint, :, :, :) = ftSlice;
%     end
% end
% 
% 
% %% -------------------------------------------------------------------------
% % getImageTrajectory
% %% -------------------------------------------------------------------------
% function [xCoordinates, yCoordinates, imageTrajectory] = getImageTrajectory(MRSIStruct)
%     xCoordinates = getCoordinates(MRSIStruct, 'x');
%     yCoordinates = getCoordinates(MRSIStruct, 'y');
% 
%     [xx, yy]      = meshgrid(xCoordinates, yCoordinates);
%     imageTrajectory = [xx(:), yy(:)];
% end
% 
% 
% %% -------------------------------------------------------------------------
% % fastFourierTransformTime: same stuff
% %% -------------------------------------------------------------------------
% function MRSIStruct = fastFourierTransformTime(MRSIStruct)
%     data = getData(MRSIStruct);
%     timeDimension = getDimension(MRSIStruct, 't');
% 
%     % Fourier transform in the spectral domain (time -> frequency)
%     data = fftshift(ifft(data, [], timeDimension), timeDimension);
%     MRSIStruct = setData(MRSIStruct, data);
% 
%     % Compute PPM axis
%     ppm = calculatePPM(MRSIStruct);
%     if strcmp(MRSIStruct.nucleus,'1H')
%         ppm = ppm + 4.65;
%     end
%     MRSIStruct = setPPM(MRSIStruct, ppm);
% 
%     % Update flags
%     MRSIStruct = setFlags(MRSIStruct, 'spectralFT', true);
% 
%     % ---- NEW: Rename dimension label from 't' to 'f' ----
%     fDim = getDimension(MRSIStruct, 't');  % Get dimension index for 't'
%     MRSIStruct = setDimension(MRSIStruct, 'f', fDim);  % Set same index for 'f'
%     MRSIStruct = setDimension(MRSIStruct, 't', 0);      % Clear old 't' label
% end
% 
% 
% %% -------------------------------------------------------------------------
% % calculatePPM: same stuff
% %% -------------------------------------------------------------------------
% function ppmVals = calculatePPM(MRSIStruct)
%     gammaVal      = MRSIStruct.gamma;
%     spectralWidth = getSpectralWidth(MRSIStruct);
%     timeSize      = getSizeFromDimensions(MRSIStruct, {'t'});
% 
%     step       = spectralWidth / timeSize;
%     lowerBound = -spectralWidth/2 + step/2;
%     upperBound =  spectralWidth/2 - step/2;
% 
%     frequencyArray = lowerBound : step : upperBound;
%     ppmVals        = -frequencyArray / (MRSIStruct.Bo * gammaVal);
% end
% 
% 
% %% -------------------------------------------------------------------------
% % calculateSpectralValues: same stuff
% %% -------------------------------------------------------------------------
% function MRSIStruct = calculateSpectralValues(MRSIStruct, kPtsPerCycle, NPtemporal)
%     spectralDwellTime = calculateSpectralDwellTime(MRSIStruct, kPtsPerCycle);
%     spectralWidth     = 1 / spectralDwellTime;
%     spectralTime      = calculateSpectralTime(spectralDwellTime, NPtemporal);
% 
%     MRSIStruct = setSpectralWidth(MRSIStruct, spectralWidth);
%     MRSIStruct = setSpectralDwellTime(MRSIStruct, spectralDwellTime);
%     MRSIStruct = setSpectralTime(MRSIStruct, spectralTime);
% end
% 
% 
% function dt = calculateSpectralDwellTime(MRSIStruct, spatialPoints)
%     adcDt = getAdcDwellTime(MRSIStruct);
%     dt    = spatialPoints * adcDt;
% end
% 
% function tAxis = calculateSpectralTime(spectralDwellTime, spatialPoints)
%     tAxis = 0:spectralDwellTime:spectralDwellTime*(spatialPoints - 1);
% end
% 
% 
% %% -------------------------------------------------------------------------
% % halfPixelShift: same stuff
% %% -------------------------------------------------------------------------
% function MRSIStruct = halfPixelShift(MRSIStruct)
%     kx = getCoordinates(MRSIStruct, 'kx');
%     ky = getCoordinates(MRSIStruct, 'ky');
% 
%     halfPixelX = getVoxSize(MRSIStruct, 'x')/2;
%     halfPixelY = getVoxSize(MRSIStruct, 'y')/2;
% 
%     kShift = kx*halfPixelX + ky'*halfPixelY;
% 
%     [MRSIStruct, prevPerm, prevSz] = reshapeDimensions(MRSIStruct, {'ky','kx'});
%     data = getData(MRSIStruct);
% 
%     data = data .* exp(-1i*2*pi*kShift);
% 
%     MRSIStruct = setData(MRSIStruct, data);
%     MRSIStruct = reshapeBack(MRSIStruct, prevPerm, prevSz);
% end

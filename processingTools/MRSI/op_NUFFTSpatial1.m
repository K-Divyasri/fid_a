function ftSpatial = op_NUFFTSpatial1(dComp, kFile_path)

    disp('=== START: op_NUFFTSpatial ===');

    %% Load and normalize k-trajectory 
    disp(['Loading k-file: ', kFile_path]);
    [kTable, ~] = readKFile(kFile_path);
    kCoords = kTable{:, 2:3};
    disp(['Original kCoords shape = ', mat2str(size(kCoords))]);

    kCoords = (kCoords / max(abs(kCoords(:)))) * pi; % normalize [-pi, pi] dont forget to normalize!!!
    disp(['Normalized kCoords shape = ', mat2str(size(kCoords))]);

    % AXIS SWAP FIX!! LOOK AT THISSS ERRROR #1
    % Swap columns to match [omega_y, omega_x] for Fessler's NUFFT
    kCoords = [kCoords(:,2), kCoords(:,1)];    % Now [Ky, Kx]
    disp('Swapped kCoords columns to [Ky, Kx] for [Ny, Nx] image output.');

    % Get dimensions from new data structure
    nKshot = dComp.sz(dComp.dims.kshot);
    nKpts = dComp.sz(dComp.dims.kpts);
    totalKPoints = nKshot * nKpts;
    
    disp(['nKshot = ', num2str(nKshot)]);
    disp(['nKpts = ', num2str(nKpts)]);
    disp(['totalKPoints = ', num2str(totalKPoints)]);

    % Verify k-space file matches data dimensions
    if size(kCoords, 1) ~= totalKPoints
        error('K-space file length (%d) does not match data k-space dimensions (%d)', ...
              size(kCoords, 1), totalKPoints);
    end

    %% NUFFT setup
    Nx = length(dComp.coordinates.x);
    Ny = length(dComp.coordinates.y);
    disp(['Nx = ', num2str(Nx), ', Ny = ', num2str(Ny)]);

    Nd = [Ny, Nx];    % NOTE THE FLIP: [Ny, Nx] THIS TOO!! ERROR #2 not an erro
    Jd = [6, 6]; 
    Kd = [2*Ny, 2*Nx]; 
    n_shift = Nd/2;

    disp('Initializing NUFFT...');
    st = nufft_init(kCoords, Nd, Jd, Kd, n_shift); % NUFFT init
    disp('NUFFT initialized.');

    %% Handle 5D or 4D data
    dataSize = size(dComp.data);
    disp(['dComp.data shape = ', mat2str(dataSize)]);
    nDims = ndims(dComp.data);
    
    hasAverages = (dComp.dims.averages > 0);

    if hasAverages && nDims == 5
        disp('=== 5D data branch (with averages) ===');
        % Data structure: [t, coils, averages, kshot, kpts]
        % Permute to [kshot, kpts, t, coils, averages]
        dataPermuted = permute(dComp.data, [4, 5, 1, 2, 3]);
        disp(['dataPermuted shape = ', mat2str(size(dataPermuted))]);
        
        % Reshape to [totalKPoints, t*coils*averages]
        dataFlat = reshape(dataPermuted, [totalKPoints, numel(dComp.data)/totalKPoints]);
        disp(['dataFlat shape = ', mat2str(size(dataFlat))]);

        imgFlat = nufft_adj(dataFlat, st);   % imgFlat is [Ny, Nx, t*coils*averages]
        disp(['imgFlat shape = ', mat2str(size(imgFlat))]);

        % Reshape to [Ny, Nx, t, coils, averages]
        img = reshape(imgFlat, [Ny, Nx, dComp.sz(1), dComp.sz(2), dComp.sz(3)]);

        disp(['img shape = ', mat2str(size(img))]);
        % Permute to [t, coils, averages, Ny, Nx] (if your pipeline expects [t, coils, averages, x, y], swap 4 and 5 here)
        imgData = permute(img, [3, 4, 5, 1, 2]); % [t, coils, averages, Ny, Nx]

        % Update dimensions
        dComp = setDimension(dComp, 'x', 5);
        dComp = setDimension(dComp, 'y', 4);
        dComp = setDimension(dComp, 'kshot', 0);
        dComp = setDimension(dComp, 'kpts', 0);

    elseif ~hasAverages && nDims == 4
        disp('=== 4D data branch (without averages) ===');
        % Data structure: [t, coils, kshot, kpts]
        % Permute to [kshot, kpts, t, coils]
        dataPermuted = permute(dComp.data, [3, 4, 1, 2]);
        disp(['dataPermuted shape = ', mat2str(size(dataPermuted))]);
        
        % Reshape to [totalKPoints, t*coils]
        dataFlat = reshape(dataPermuted, [totalKPoints, numel(dComp.data)/totalKPoints]);
        disp(['dataFlat shape = ', mat2str(size(dataFlat))]);

        imgFlat = nufft_adj(dataFlat, st);

        disp(['imgFlat shape = ', mat2str(size(imgFlat))]);

        % Reshape to [Ny, Nx, t, coils]
        img = reshape(imgFlat, [Ny, Nx, dComp.sz(1), dComp.sz(2)]);

        disp(['img shape = ', mat2str(size(img))]);
        % Permute to [t, coils, Ny, Nx]
        imgData = permute(img, [3, 4, 1, 2]);

        % Update dimensions
        dComp = setDimension(dComp, 'x', 4);
        dComp = setDimension(dComp, 'y', 3);
        dComp = setDimension(dComp, 'kshot', 0);
        dComp = setDimension(dComp, 'kpts', 0);

    else
        error('Unsupported data dimension: %d or averages flag inconsistent', nDims);
    end

    % Update data and flags
    dComp.data = imgData;
    dComp.sz = size(imgData);
    dComp.flags.spatialFT = 1;
    
    % Clear k-space dimensions
    dComp = setDimension(dComp, 'kx', 0);
    dComp = setDimension(dComp, 'ky', 0);
    
    % Spectral updates
    kPtsPerCycle = getKPtsPerCycle(kTable);
    disp(['Debug: kPtsPerCycle = ', mat2str(kPtsPerCycle)]);
    NPtemporal = getTemporalPts(kTable, dComp);
    disp(['Debug: NPtemporal = ', mat2str(NPtemporal)]);
    dComp = calculateSpectralValues(dComp, kPtsPerCycle, NPtemporal);
    
    ftSpatial = dComp;

    disp('=== END: op_NUFFTSpatial ===');
end
function MRSIStruct = calculateSpectralValues(MRSIStruct, kPtsPerCycle, NPtemporal)
    spectralDwellTime = calculateSpectralDwellTime(MRSIStruct, kPtsPerCycle);
    spectralWidth = 1/spectralDwellTime;
    spectralTime = calculateSpectralTime(spectralDwellTime, NPtemporal);

    MRSIStruct = setSpectralWidth(MRSIStruct, spectralWidth);
    MRSIStruct = setSpectralDwellTime(MRSIStruct, spectralDwellTime);
    MRSIStruct = setSpectralTime(MRSIStruct, spectralTime);
end

function spectralDwellTime = calculateSpectralDwellTime(MRSIStruct, spatialPoints)
    adcDwellTime = getAdcDwellTime(MRSIStruct);
    spectralDwellTime = spatialPoints * adcDwellTime;
end

function spectralTime = calculateSpectralTime(spectralDwellTime, spatialPoints)
    spectralTime = 0:spectralDwellTime:spectralDwellTime*(spatialPoints - 1);
end

function MRSIStruct = halfPixelShift(MRSIStruct)
    kx = getCoordinates(MRSIStruct, 'kx');
    ky = getCoordinates(MRSIStruct, 'ky');
    halfPixelX = getVoxSize(MRSIStruct, 'x')/2;
    halfPixelY = getVoxSize(MRSIStruct, 'y')/2;
    kShift = kx*halfPixelX + ky'*halfPixelY;

    [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'ky', 'kx'});
    data = getData(MRSIStruct);
    data = data .* exp(-1i*2*pi*kShift);
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);
end

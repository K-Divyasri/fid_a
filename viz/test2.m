final_location = 'C:\Users\divya\Downloads\linux-send\fida3 - Copy (2)\new\lcm_op\lcm_op_rosette1';
figurefoldername = 'Figures_lcm_op_rosette1';

[map, crlb, LW, SNR] = op_CSILCModelMaps(40,40,final_location,'figure_folder_name',figurefoldername);

% Visualise (example)
names = {'CrPCr','GPCPCh','NAANAAG','GluGln','Ins'};
for i=1:numel(names)
    figure('Name',names{i});
    imagesc(map.(names{i})); axis image; colorbar; colormap hot;
    title(['LCModel ' names{i}]);
end
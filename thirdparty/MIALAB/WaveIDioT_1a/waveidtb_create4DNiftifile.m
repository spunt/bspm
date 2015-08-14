function waveidttb_create4DNiftifile(outputFileName, data_fft, mat)
% Function to create 4D file

% mat = [-1.1420 0 0.1229  107.6518; ...
%     -0.0176    1.1255   -0.2803  -35.9283; ...
%     0.0922    0.2148    1.4685  -92.6349; ...
%     0         0         0    1.0000];

%  VM = spm_vol('/export/apps/linux-x86/matlab/toolboxes/spm5/templates/EPI.nii');
%  mat = VM(1).mat;


V = struct('fname',  outputFileName,...
    'dim',    [size(data_fft, 1), size(data_fft, 2),size(data_fft, 3)],...
    'dt',     [spm_type('int16'), spm_platform('bigend')],...
    'mat',    mat,...
    'pinfo',  [1 0 0]',...
    'descrip','spm:');

mx   = -Inf;
mn   = Inf;
for i=1:size(data_fft, 4)
    %dat      = V(i).private.dat(:,:,:,ind(i,1),ind(i,2));
    dat = squeeze(data_fft(:, :, :, i));
    dat      = dat(isfinite(dat));
    mx       = max(mx,max(dat(:)));
    mn       = min(mn,min(dat(:)));
end;

sf         = max(mx,-mn)/32767;
ni         = nifti;
ni.dat     = file_array(outputFileName, [V.dim(1:3), size(data_fft, 4)], 'INT16-BE',0,sf,0);
ni.mat     = V.mat; %N(1).mat;
ni.mat0    = V.mat; %N(1).mat;
ni.descrip = V.descrip;
disp(['Writing 4D file: ', outputFileName]);
create(ni);
for i=1:size(ni.dat,4),
    ni.dat(:,:,:,i) = data_fft(:, :, :, i); %N(i).dat(:,:,:,ind(i,1),ind(i,2));
    spm_get_space([ni.dat.fname ',' num2str(i)], V.mat);
%     i
end;
disp('Done writing 4D file');
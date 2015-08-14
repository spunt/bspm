function [mask ,x , y , z , w ] = waveidtb_generateMask(ParameterInfo)

% Find the number of timepoints (if NII) else use the one entered by user (IF
% img)
if strcmpi(ParameterInfo.fileFormat,'nii')&& strcmpi(ParameterInfo.fileNames{1}(1,end-4:end-2),'img')
    V = spm_vol(ParameterInfo.fileNames{1});
    firstVols =  unique(strtok(ParameterInfo.fileNames,'.'));
    % Create the filenames for first volumes only
    fname = cell2mat([firstVols repmat({'.nii,1'},length(firstVols),1)]);
    
elseif strcmpi(ParameterInfo.fileFormat,'nii')
    V = spm_vol(ParameterInfo.fileNames{1});
    firstVols =  unique(strtok(ParameterInfo.fileNames,'.'));
    % Create the filenames for first volumes only
    fname = cell2mat([firstVols repmat({'.nii,1'},length(firstVols),1)]);
    
else
    V = spm_vol(ParameterInfo.fileNames{1});
    firstVols =  unique(strtok(ParameterInfo.fileNames,'.'));
    fname = cell2mat([firstVols repmat({['.',ParameterInfo.fileFormat ',1']},length(firstVols),1)]);
end;

%%%%%%%%%%%%%%%% READING THE DATA (Can use GIFT's method) %%%%%%%%%%%%
fprintf('\nGenerating mask using first volume from each Subject/Session.\n');
img = spm_read_vols( spm_vol( fname ) );
[x y z,w] = size( img ); % redefine size

mask = ones(x,y,z);

for i = 1:size(img,4)
    %% Estimate the size of the
    temp = squeeze(img(:,:,:,i));

    %% Calculate the Mask for the First subject and use that for all
    % subjects
      
    threshold= 0.80*mean( temp ( ~isnan(temp) ) );
    maskMean = double(abs( temp >= ( threshold .* ones( size( temp ) ) ) ));
    
    % Fill the isolated holes in the intracrennial voxels
    for i = 1:size(maskMean,3);
        maskMean(:,:,i) = double( bwfill( bwmorph(maskMean(:,:,i)>0,'clean') ,'holes'));
    end;
    mask = and(maskMean, mask);
    clear maskMean temp;
end;

if length(V.private.dat.dim) >= 4
    w = V.private.dat.dim(4);
elseif w>=1
    fprintf('Done.\n');
    return;
else
    w = 1;
end;



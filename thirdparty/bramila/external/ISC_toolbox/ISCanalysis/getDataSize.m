function [siz,flag] = getDataSize(fileName,fileFormat)

flag = true;

if strcmp(fileFormat,'nii') % load nii-file
    I = load_nii_hdr(fileName);
    I = I.dime.dim;%uint8(I.img);
    if I(1) ~= 4
        if nargin == 1
            error('fMRI data must be 4D!')            
        else
           flag = false;
           siz = NaN*ones(1,4);
           return
        end
    end
    siz = I(2:5);
elseif strcmp(fileFormat,'mat') % load mat-file
    q = whos('-file',fileName);
    siz = q(1).size;
    if length(siz) ~= 4
        if nargin == 1
            error('fMRI data must be 4D!')
        else
            flag = false;
            siz = NaN*ones(1,4);
            return
        end
    end
else
    error('Preprocessed data must be of type .nii or .mat!!')
end

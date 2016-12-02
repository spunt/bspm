function bnii_update_dicm2nii
% function SUCCESS = bnii_update_dicm2nii
%
    baseurl     = 'http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions';
    url         = 'http://www.mathworks.com/matlabcentral/fileexchange/42997-dicom-to-nifti-converter--nifti-tool-and-viewer';
    html        = urlread(url);
    str         = regexp(html, '/\d{5,5}/versions/\d+/', 'match');
    dlurl       = sprintf('%s%sdownload/zip', baseurl, str{1});
    outdir      = fileparts(which('dicm2nii.m'));
    locpath     = fullfile(outdir, 'dicm2nii.zip');
    fprintf('\nDownloading from %s\nWriting to %s', dlurl, locpath);
    [F, STATUS] = urlwrite(dlurl, locpath);
    if STATUS
        fprintf('\nDownload complete. Unzipping package...', dlurl, locpath);
        unzip(locpath, outdir);
        fprintf(' FINISHED\n');
    else
        fprintf('\nDownload failed!\n');
    end

end


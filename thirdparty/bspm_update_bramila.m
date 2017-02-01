function bspm_update_bramila
% BSPM_UPDATE_BRAMILA
%
dlurl = 'https://git.becs.aalto.fi/bml/bramila/repository/archive.zip';
outdir      = fileparts(which('bramila_dvars.m'));
locpath     = fullfile(outdir, 'archive.zip');
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


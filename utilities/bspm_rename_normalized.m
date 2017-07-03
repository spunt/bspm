function bspm_rename_normalized(in)
    % bspm_rename_normalized(in)
    if nargin < 1, mfile_showhelp; return; end
    if ischar(in), in = cellstr(in); end
    [pth, fname, fext] = cellfun(@fileparts, in, 'Unif', false);
    hdr = spm_vol(char(in));
    des = {hdr.descrip}';
    nimg = length(in);
    out = in;
    for i = 1:nimg
        [p,n,e] = fileparts(in{i});
        h = spm_vol(in{i});
        res = abs(h(1).mat(1));
        fwhm = regexp(h(1).descrip, '\dx', 'once', 'match');
        fwhm = fwhm(1);
        out{i} = fullfile(p, strcat(reegexprep(n, '^sw', sprintf('s%sw%d',fwhm, res)), e));
    end
    cellfun(@movefile, in, out);
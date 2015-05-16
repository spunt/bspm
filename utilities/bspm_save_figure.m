function im = bspm_save_figure(name, varargin)
% BOB_SAVE_FIGURE
% 
% USAGE: im = bspm_save_figure(name, varargin);
%
% ARGUMENTS
%   name: output name
% VARARGIN
%   1 - fmt: [jpg], png, pdf, eps, tif, bmp
%   2 - renderer: [opengl], painters
%   2 - crop: [0]=nocrop, 1=crop margins
%   3 - magnify: factor to magnify figure size by (integer)
%   4 - transparent: [0]=no, 1=make backg transparent (png, pdf, eps only)
%   5 - ppi: pixels per inch, or resolution (e.g., 300)
%
% =================================================================================
def = {     'fmt',          'jpg',      ...
            'renderer',     'opengl',   ...
            'crop',         0,          ...
            'magnify',      1,          ...
            'ppi',          100,        ...
            'transparent',  0,          ...
            'timestamp',    1           ...
            }; 
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
fmt = strcat('-', fmt);
if iscell(name), name = char(name); end
[fp, fn, fe] = fileparts(name); 
if isempty(fe), fe = ['.' fmt]; end
if isempty(fp), fp = pwd; end
if timestamp, fn = [fn '_' sprintf('%s_%s', datestr(now,'mm-DD-YYYY'), lower(strtrim(datestr(now, 'HHMMpm'))))]; end
name = fullfile(fp, [fn fe]);
mgn = sprintf('-m%d', magnify); 
rnd = strcat('-', renderer);
res = sprintf('-r%d', ppi); 
opt = {fmt mgn rnd res};
if transparent, opt{end+1} = '-transparent'; end
if ~crop, opt{end+1} = '-nocrop'; end
if nargout==1
    im = export_fig(name, opt{:}); 
else
    export_fig(name, opt{:}); 
end
end
 
 
 
 

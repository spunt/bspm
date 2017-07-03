function matlabbatch = bspm_contrasts(analysis_dir, weights, delete_tag, repl_tag, connames)
% BSPM_CONTRASTS
%
%   USAGE: bspm_contrasts(analysis_dir, weights, delete_tag, repl_tag, connames)
%
%   ARGUMENTS:
%       analysis_dir: directory containing SPM.mat
%       weights: contrast weights
%       delete_tag: tag for deleting existing contrasts: 0 for no (default), 1 for yes
%       repl_tag: tag for replicating across sessions: 0 for no, 1 for yes (default)
%       connames: custom contrast names
%

% -------------------------------- Copyright (C) 2014 --------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<2, mfile_showhelp; return; end
if nargin<5, connames = []; end
if nargin<4, repl_tag = 1; end
if nargin<3, delete_tag = 0; end
if ischar(analysis_dir), analysis_dir = cellstr(analysis_dir); end
if repl_tag, repl_choice = 'repl'; else repl_choice = 'none'; end
for s = 1:length(analysis_dir)
    
spmmat = [analysis_dir{s} filesep 'SPM.mat'];
ncon = size(weights,1);

if isempty(connames)
    
    tmp = load(spmmat);
    regnames = tmp.SPM.xX.name;
    % clean up names
    for r = 1:length(tmp.SPM.Sess)
        string = sprintf('Sn\\(%d\\) ', r);
        regnames = regexprep(regnames,string,'');
    end
    for r = 1:length(regnames)
        n = regnames{r};
        n(strfind(n,'*'):end) = [];
        n(strfind(n,'^'):end) = [];
        regnames{r} = n;
    end
    condnames = regnames;

    % build contrast names
    
    for c = 1:ncon
        w = weights(c,:);
        % build name
        posidx = [];
        posidx = find(w>0);
        negidx = [];
        negidx = find(w<0);
        if isempty(negidx)
            name = strcat(condnames{posidx});
        elseif isempty(posidx)
            name = ['INV_' strcat(condnames{negidx})];
        else
            name = [strcat(condnames{posidx}) '_-_' strcat(condnames{negidx})];
        end
        connames{c} = name;
    end
    
else
    
    if ischar(connames), connames = cellstr(connames); end
    if length(connames)~=size(weights, 1)
        disp('Number of contrast names does not equal number of weight vectors. Exiting...');
        matlabbatch = [];
        return;
    end
    
end

% build job
matlabbatch{s}.spm.stats.con.spmmat{1} = spmmat;
matlabbatch{s}.spm.stats.con.delete = delete_tag;
for c = 1:ncon
    matlabbatch{s}.spm.stats.con.consess{c}.tcon.name       = connames{c};
    matlabbatch{s}.spm.stats.con.consess{c}.tcon.weights    = weights(c,:); % SPM12
    matlabbatch{s}.spm.stats.con.consess{c}.tcon.convec     = weights(c,:); % SPM8
    matlabbatch{s}.spm.stats.con.consess{c}.tcon.sessrep    = repl_choice;
end

end

% run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end

end

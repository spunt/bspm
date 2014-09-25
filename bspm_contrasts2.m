function bspm_contrasts2(input)
% BSPM_CONTRASTS
%
%   USAGE: bspm_contrasts(analysis_dir, weights, delete_tag, repl_tag)
%
%   ARGUMENTS:
%       analysis_dir: directory containing SPM.mat
%       weights: contrast weights
%       delete_tag: tag for deleting existing contrasts: 0 for no (default), 1 for yes
%       repl_tag: tag for replicating across sessions: 0 for no, 1 for yes (default)
%

% -------------------------------- Copyright (C) 2014 --------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin < 1, error('No input!'); end
if ~isfield(input, 'delete_tag'), input.delete_tag = 0; end
if ~isfield(input, 'repl_tag'), input.repl_tag = 1; end
fn = {'analysis_dir' 'weights' 'delete_tag' 'repl_tag'};
[status, msg] = checkfields(input, fn);
if ~status, error(msg); end
analysis_dir = input.analysis_dir; 
weights = input.weights; 
repl_tag = input.repl_tag; 
delete_tag = input.delete_tag; 
if ischar(analysis_dir), analysis_dir = cellstr(analysis_dir); end
if repl_tag, repl_choice = 'repl'; else repl_choice = 'none'; end

for s = 1:length(analysis_dir)

spmmat = [analysis_dir{s} filesep 'SPM.mat'];
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
ncon = size(weights,1);
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

% build job
matlabbatch{1}.spm.stats.con.spmmat{1} = spmmat;
matlabbatch{1}.spm.stats.con.delete = delete_tag;   
for c = 1:ncon
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = connames{c};
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = weights(c,:);
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = repl_choice;
end

% run job
spm('defaults','fmri'); spm_jobman('initcfg');         
spm_jobman('run',matlabbatch);

end

end

 
 
 
 

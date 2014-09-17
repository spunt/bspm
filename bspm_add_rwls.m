function bspm_add_rwls(nuisancefile, goodspmmat)
% BSPM_ADD_RWLS
%
%   'USAGE = bspm_add_rwls(nuisancefile, goodspmmat)'
%
%       nuisance_file:      txt file with nuisance regressors 
%       goodspmmat:         SPM.mat with rWLS regressor(s)
%

% --------------------------------- Copyright (C) 2014 ---------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<2, error('USAGE = bspm_add_rwls(nuisancefile, goodspmmat)'); end
if ischar(nuisancefile) nuisancefile = cellstr(nuisancefile); end

%% GET RWLS
[x,xn,sidx]  = bspm_get_design(goodspmmat, 1);
nsess   = length(unique(sidx));
if length(nuisancefile)~=nsess, error('Too few nuisance files provided!'); end
rwls    = x(:,strcmp(xn, 'rWLS'));

%% LOOP OVER NUISANCE FILES
for i = 1:nsess
    
    orig    = load(nuisancefile{i});
    new     = [orig rwls(sidx==i, i)];
    save(nuisancefile{i}, 'new', '-ascii'); 
    
end

 
 
 
 

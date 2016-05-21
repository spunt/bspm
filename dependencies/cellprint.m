function outcell = cellprint(incell, rmpath, forcenumber, nonumber)
% CELLPRINT Show cell array
%
%       USAGE: cellprint(incell, [rmpath], [forcenumber], [nonumber])
%           
% -------------------------------------------------------------------
if nargin<4, nonumber = 0; end
if nargin<3, forcenumber = 0; end
if nargin<2, rmpath = 1; end
if nargin<1, mfile_showhelp; return; end
if ischar(incell), incell = cellstr(incell); end
if rmpath, [~,incell,~] = cellfun(@fileparts, incell, 'Unif', false); end
if size(incell,1)==1, incell = incell'; end
[nrow, ncol] = size(incell); 
if ~nonumber & (ncol==1 | forcenumber)
    % | Number
    N           = cell(length(incell), 1);
    N(:,1)      = num2cell(1:length(incell));
    str         = sprintf('%%0%dd', max(cellfun('length', cellfun(@num2str, N, 'Unif', false))));
    for i = 1:length(incell), N{i} = sprintf(str, i); end
    incell(:,1) = strcat({' - '}, incell(:,1));
    incell(:,1) = strcat(cellstr(N), incell(:,1));
end
if ncol>1
    ln      = cellfun('length', incell);
    npad    = repmat(max(ln), nrow, 1) - ln + 1;
    for i = 1:nrow
        for c = 1:ncol
            incell{i,c} = [incell{i,c} repmat(' ', 1, npad(i,c))]; 
        end
    end
end
for i = 1:nrow, fprintf(['\n ' repmat('%s', 1, ncol)], incell{i,:}); end
fprintf('\n\n');
if nargin > 0, outcell = incell; end

end
function ind = ipctb_ica_good_cells( mycell )

if ~iscell(mycell)
    mycell = {mycell};
end

for j=1:length(mycell)
    ind(j)=~isempty(mycell{j});
end

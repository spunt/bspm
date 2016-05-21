function hprop = arrayget(harray, propname) 
% ARRAYGET Get handles for multiple axes
%
% USAGE: h = arrayget(array, propname)
%
% ==============================================
if nargin<2, error('USAGE: h = arrayget(array, propname)'); end
if size(harray, 1)==1, harray  = harray'; end
hprop = get(harray, propname); 
hprop = [hprop{:}]; 
% hprop = cell2mat(arrayfun(@get, harray, repmat({propname}, length(harray), 1)));
end
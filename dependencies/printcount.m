function thecount = printcount(i, n)
% PRINTCOUNT Print progress count to command window
%
%  USAGE: thecount = printcount(i, n)
% __________________________________________________________________________
%  INPUTS
% 	i: current count   
% 	n: total count
% __________________________________________________________________________
%  EXAMPLE USAGE
%   myarray = 1:42; 
%   fprintf('Count: '); 
%   for i = 1:length(myarray) 
%       printcount(i, length(myarray));
%       pause(.1);
%   end
%   fprintf('\n'); 
%       

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-03-17
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 2, disp('USAGE: thecount = printcount(i, n)'); return; end
nd  = length(num2str(n)); 
f   = sprintf('%%0%dd', nd); 
if nargout, thecount = sprintf([f '/' f], i, n); return; end
if i>1
    fprintf([repmat('\b', 1, nd*2+1) f '/' f], i, n);
else
    fprintf([f '/' f], i, n);
end
if i==n, fprintf('\n'); end


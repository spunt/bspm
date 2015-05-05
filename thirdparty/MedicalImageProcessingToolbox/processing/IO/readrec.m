function data = readrec(fname,par,format,data_type,n0)

% y = read_rec(FNAME,PAR); 
% y = read_rec(FNAME,PAR,FORMAT);
% y = read_rec(FNAME,PAR,FORMAT,NPH);  NPH : phases to read
%
% Read .REC file FNAME using specified PAR info and FORMAT (default uint16)
% If parameters from the PAR file are given, data is reshaped to 
% [nx ny slices phases dynamics echoes image_type]

if (nargin < 3), 
    format = 'uint16';
    %format = 'uint32';
end
if (nargin == 1),
  par = []; 
end

pid = fopen(fname,'r');
if (pid == -1),
    error('could not open file...');
end
data = fread(pid,format);
%data = single(data);
fclose(pid);

if ~isempty(par),
  % reshape the data [nx ny slices phases dynamics echoes image_type]
  nx = par.dim(1); ny = par.dim(2); nz = par.slice; 
  phases = par.phases; dyn = max(par.dyn);

  types =unique(par.ty);
  %imtype = types(end); % type = 0 for modulus image
  %imtype = numel(types); % type = 0 for modulus image
  imtype = 1; % number of images that this data file contains
  %imtype = 2; % number of images that this data file contains
  %imtype = 3; % number of images that this data file contains
  
  echoes = par.echoes;
  
  if imtype==0
    NZ = length(data)/nx/ny/phases/dyn/echoes;
    data = reshape(data,[nx ny nz phases dyn echoes ]);
  else
    NZ = length(data)/nx/ny/phases/dyn/echoes/imtype;
 
      if (NZ > nz),
        fprintf('kz oversampling factor: %1.4f',NZ/nz);
        nz = NZ;
      end
     % data = permute(reshape(data,[ nx ny  phases  nz dyn echoes imtype]),[1 2 4 3 5 6 7]);
      
      data=reshape(data,[ nx ny   nz phases dyn echoes imtype]);
  end
  
  %data = reshape(data,[ny nx nz phases dyn echoes imtype]);

  data = single(data);
  
%    % scale data using PAR info
%   idx = sub2ind([nz phases dyn echoes imtype],par.sl,par.ph,par.dyn,par.ec,par.ty);
%   par.ri(idx) = single(par.ri); par.ri = reshape(par.ri,[1 1 nz phases dyn echoes imtype]);
%   par.rs(idx) = single(par.rs); par.rs = reshape(par.rs,[1 1 nz phases dyn echoes imtype]);
%   par.ss(idx) = single(par.ss); par.ss = reshape(par.ss,[1 1 nz phases dyn echoes imtype]);
%   data = (data .* repmat(par.rs,[nx ny 1 1 1 1 1]) + repmat(par.ri,[nx ny 1 1 1 1 1])) ./ repmat((par.rs .* par.ss),[nx ny 1 1 1 1 1 1]);
  
end


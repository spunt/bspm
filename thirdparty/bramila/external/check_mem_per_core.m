function memreq_line = check_mem_per_core(inputfile,mb_per_TR)
% Generic function to handle variable memory requirements.
% If job memory usage depends significantly on file size, it may be good to
% optimize it. Smaller memory requirements make jobs get their queue much
% faster. Use this as template.
%
% Input
% -inputfile : path to 4D image file, whose size will be checked
% -estimates : 
% Output
% -memreq_line : line for SLURM job file that was formatted accordingly
%
if nargin == 1
    % Case of preprocessor,from experience, preprocessing requires around 35gb memory if the file is 1055 TR
    % Adjust a little bit just in case
	% Enrico increased:
    mb_per_TR = 60000/1055;
end
hdr = load_nii_hdr(inputfile);
nvol = hdr.dime.dim(5);
% [~,nvol]=system(sprintf('fslnvols %s',inputfile)); nvol = str2double(nvol);
memreq = ceil(nvol * mb_per_TR);
memreq_line = sprintf('#SBATCH --mem-per-cpu=%i',memreq);

function [V, sizes, origin,scales] =read_gipl2(filename)
% function for reading header of  Guys Image Processing Lab (Gipl) volume file
%
% [V, sizes, origin,scales] = gipl_read_header(filename);
%
%
% Copied from the package ReadData3D_version1h


fid=fopen(filename,'rb','ieee-be');
%fid=fopen(filename,'r','native');
if(fid<0)
    fprintf('could not open file %s\n',filename);
    return
end

trans_type{1}='binary'; trans_type{7}='char'; trans_type{8}='uchar'; 
trans_type{15}='short'; trans_type{16}='ushort'; trans_type{31}='uint'; 
trans_type{32}='int';  trans_type{64}='float';trans_type{65}='double'; 
trans_type{144}='C_short';trans_type{160}='C_int';  trans_type{192}='C_float'; 
trans_type{193}='C_double'; trans_type{200}='surface'; trans_type{201}='polygon';

trans_orien{0+1}='UNDEFINED'; trans_orien{1+1}='UNDEFINED_PROJECTION'; 
trans_orien{2+1}='AP_PROJECTION';  trans_orien{3+1}='LATERAL_PROJECTION'; 
trans_orien{4+1}='OBLIQUE_PROJECTION'; trans_orien{8+1}='UNDEFINED_TOMO'; 
trans_orien{9+1}='AXIAL'; trans_orien{10+1}='CORONAL'; 
trans_orien{11+1}='SAGITTAL'; trans_orien{12+1}='OBLIQUE_TOMO';

offset=256; % header size

%get the file size
fseek(fid,0,'eof');
fsize = ftell(fid); 
fseek(fid,0,'bof');

sizes=fread(fid,4,'ushort')';
if(sizes(4)==1), maxdim=3; else maxdim=4; end
sizes=sizes(1:maxdim);
image_type=fread(fid,1,'ushort');
scales=fread(fid,4,'float')';
scales=scales(1:maxdim);
patient=fread(fid,80, 'uint8=>char')';
matrix=fread(fid,20,'float')';
orientation=fread(fid,1, 'uint8')';
par2=fread(fid,1, 'uint8')';
voxmin=fread(fid,1,'double');
voxmax=fread(fid,1,'double');
origin=fread(fid,4,'double')';
origin=origin(1:maxdim);
pixval_offset=fread(fid,1,'float');
pixval_cal=fread(fid,1,'float');
interslicegap=fread(fid,1,'float');
user_def2=fread(fid,1,'float');
magic_number= fread(fid,1,'uint32');
%if (magic_number~=4026526128), error('file corrupt - or not big endian'); end
fclose('all');

info.Filename=filename;
info.FileSize=fsize;
info.Dimensions=sizes;
info.PixelDimensions=scales;
info.ImageType=image_type;
info.Patient=patient;
info.Matrix=matrix;
info.Orientation=orientation;
info.VoxelMin=voxmin;
info.VoxelMax=voxmax;
info.Origing=origin;
info.PixValOffset=pixval_offset;
info.PixValCal=pixval_cal;
info.InterSliceGap=interslicegap;
info.UserDef2=user_def2;
info.Par2=par2;
info.Offset=offset;

V = gipl_read_volume(info);


end

function V = gipl_read_volume(info)
% function for reading volume of Guys Image Processing Lab (Gipl) volume file
% 
% volume = gipl_read_volume(file-header)
%
% examples:
% 1: info = gipl_read_header()
%    V = gipl_read_volume(info);
%    imshow(squeeze(V(:,:,round(end/2))),[]);
%
% 2: V = gipl_read_volume('test.gipl');

if(~isstruct(info)) info=gipl_read_header(info); end

% Open gipl file
fid=fopen(info.Filename','rb','ieee-be');

  % Seek volume data start
  if(info.ImageType==1), voxelbits=1; end
  if(info.ImageType==7||info.ImageType==8), voxelbits=8; end
  if(info.ImageType==15||info.ImageType==16), voxelbits=16; end
  if(info.ImageType==31||info.ImageType==32||info.ImageType==64), voxelbits=32; end
  if(info.ImageType==65), voxelbits=64; end
  
  datasize=prod(info.Dimensions)*(voxelbits/8);
  fsize=info.FileSize;
  fseek(fid,fsize-datasize,'bof');

  % Read Volume data
  volsize(1:3)=info.Dimensions;

  if(info.ImageType==1), V = logical(fread(fid,datasize,'bit1')); end
  if(info.ImageType==7), V = int8(fread(fid,datasize,'char')); end
  if(info.ImageType==8), V = uint8(fread(fid,datasize,'uchar')); end
  if(info.ImageType==15), V = int16(fread(fid,datasize,'short')); end 
  if(info.ImageType==16), V = uint16(fread(fid,datasize,'ushort')); end
  if(info.ImageType==31), V = uint32(fread(fid,datasize,'uint')); end
  if(info.ImageType==32), V = int32(fread(fid,datasize,'int')); end
  if(info.ImageType==64), V = single(fread(fid,datasize,'float')); end 
  if(info.ImageType==65), V = double(fread(fid,datasize,'double')); end 

fclose(fid);

% Reshape the volume data to the right dimensions
V = reshape(V,volsize);

end



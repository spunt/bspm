function out = setOrientation(im,mat,varargin)
% out = setOrientation(im,mat)
% out = resampleImage(im,mat, options)
%
%   Options:
%       'invert',use inverse matrix
%
%   Works for ImageType and VectorImageType



invert=false;
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'invert'))
        invert=true;
    end
    
end
%----------------------------

imageMatrix = eye(4);
imageMatrix(1:3,1:3) = im.orientation(1:3,1:3);
imageMatrix(1:3,4) = im.origin(1:3);

if ~invert
    matrix_use = inv(mat);
else
     matrix_use = mat;
end

orientation_matrix = matrix_use * imageMatrix;

newOrigin = matrix_use*[im.origin(1:3) ; 1];
no = im.origin;
no(1:3)=newOrigin(1:3);

om = im.orientation;
om(1:3,1:3) = orientation_matrix(1:3,1:3);

if (isa(im,'VectorImageType') || isfield(im,'datax'))
    out = VectorImageType(im.size,no,im.spacing,om);
    vectorsIn = [im.datax(:) im.datay(:) im.dataz(:)]';
    
    vectorsOut = matrix_use(1:3,1:3)*vectorsIn ;
    
    out.datax(:) = vectorsOut(1,:)';
    out.datay(:) = vectorsOut(2,:)';
    out.dataz(:) = vectorsOut(3,:)';
    
    % TOSO transform vectors also
else
    out = ImageType(im.size,no,im.spacing,om);
    out.data = im.data;
end





end
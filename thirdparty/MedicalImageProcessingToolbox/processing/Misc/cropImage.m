function cropped_image=cropImage(imag,boundaries)
%CROPIMAGE crops an image given a bounding box in world coordinates
% cropped=cropImage(imag,boundaries)
%
%
% boundaries is a column vector vector
    
    epsilon=1E-03;
    boundaries_input = imag.GetBounds();
    nd = numel(boundaries_input)/2;
    % convert  boundaries (8 corners) into voxel index of the input image
    
    npoints = 2^nd;
    
    
    str1 = ''; str2 = ''; str3='';
    for idx=1:nd
        str1=[str1 'i' num2str(idx) ' '];
        str2=[str2 '[' num2str(2*(idx-1)+1) ':' num2str(2*(idx-1)+2) '],'];
        str3=[str3 'i' num2str(idx) '(:) '];
    end
    eval([ '[' str1(1:end-1) ']=ndgrid(' str2(1:end-1) '); indices=[' str3(1:end-1) ']; clear ' str1 ';']);

    points = boundaries(indices); 
        
    continuousIndices = imag.GetContinuousIndex(points');
    
    index1 = max([floor(min(continuousIndices,[],2)+epsilon) ones(nd,1)],[],2);
    index2 = min([ceil(max(continuousIndices,[],2 )-epsilon)  imag.size ],[],2);
    
    % If the output extent is larger than the input extent,more voxels have to be added    
    
  str1 = ''; 
    for idx=1:nd
        str1=[str1 'max(index1(' num2str(idx) '),1):min(index2(' num2str(idx) '),size(imag.data,' num2str(idx) ')),' ];
    end
    eval([ 'cropped.data = imag.data(' str1(1:end-1) ');']);
   
    cropped.spacing = imag.spacing;
    %origin must be changed now
    %change this origin  for the first true point after this one!
    
    cropped.origin = imag.GetPosition(index1);
    
    cropped_image = ImageType(size(cropped.data),cropped.origin,cropped.spacing,imag.orientation);
    cropped_image.data = cropped.data;
    
if ( numel(imag.findprop('datax'))>0) % vector image
    clear cropped_image;
    cropped_image = VectorImageType(size(cropped.data),cropped.origin,cropped.spacing,imag.orientation);

      eval([ 'cropped_image.datax = imag.datax(' str1(1:end-1) ');']);
      eval([ 'cropped_image.datay = imag.datay(' str1(1:end-1) ');']);
      eval([ 'cropped_image.dataz = imag.dataz(' str1(1:end-1) ');']);

end
    
end

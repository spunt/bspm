function outputImage = transformImage(matrix, inputImg, referenceImage,varargin)
% matrix is the transform from target to source!
%
%   'adaptSize' set to true will change output origin and size to cover all the input image
%   information
%   'adaptSizeToBoth' set to true will change output origin and size to
%   cover all the bounds of both input and output images


 mode = 'NN'; % NN
 adaptSize=false;
 adaptSizeToBoth=false;
 crop=false;
 totalBounds=[];
 MAX_CHUNK_SIZE = 50;
 for i=1:size(varargin,2)    
    if (strcmp(varargin{i},'debug'))
        debug = true;
    elseif (strcmp(varargin{i},'interpolation'))
        mode= varargin{i+1};
        i=i+1;
     elseif (strcmp(varargin{i},'maxChunkSize'))
        MAX_CHUNK_SIZE = varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'adaptSize'))
        adaptSize= true;
    elseif (strcmp(varargin{i},'adaptSizeToBoth'))
        adaptSizeToBoth= true;
    elseif (strcmp(varargin{i},'customBounds'))
        adaptSizeToBoth= true;
        totalBounds= varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'crop'))
        crop= true;
    end
    
 end
 

    if (adaptSizeToBoth)
        % compute output bounds
        if (numel(totalBounds)==0)
            boundsSource = inputImg.GetBounds();
            boundsTarget = referenceImage.GetBounds();        
            totalBounds = zeros(6,1);
            totalBounds([1 3 5])=min( [ boundsSource([1 3 5]) boundsTarget([1 3 5]) ],[],2);
            totalBounds([2 4 6])=max([ boundsSource([2 4 6]) boundsTarget([2 4 6]) ],[],2);
        end
        % compute the eight corners
        c_transformed(:,1) = matrix \ [totalBounds([1 3 5]) ;1];
        c_transformed(:,2) = matrix \ [totalBounds([1 3 6]) ;1];
        c_transformed(:,3) = matrix \ [totalBounds([1 4 5]) ;1];
        c_transformed(:,4) = matrix \ [totalBounds([1 4 6]) ;1];
        c_transformed(:,5) = matrix \ [totalBounds([2 3 5]) ;1];
        c_transformed(:,6) = matrix \ [totalBounds([2 3 6]) ;1];
        c_transformed(:,7) = matrix \ [totalBounds([2 4 5]) ;1];
        c_transformed(:,8) = matrix \ [totalBounds([2 4 6]) ;1];
        
        
        % new bounds
        m=min(c_transformed,[],2);
        M=max(c_transformed,[],2);
        if (referenceImage==referenceImage)
            newSpacing = referenceImage.spacing;
        else
            newSpacing= inputImg.spacing;
        end
        newSize = ceil((M(1:3)-m(1:3))./newSpacing)+1;
        newOrigin = m(1:3);
        newOrigin = newOrigin - rem(newOrigin-inputImg.origin,newSpacing)-newSpacing;
        
        outputImage = ImageType(newSize,newOrigin,newSpacing,eye(3));
    elseif (adaptSize)
         % compute output bounds
        boundsSource = inputImg.GetBounds();
        % compute the eight corners
        c_transformed(:,1) = matrix \ [boundsSource([1 3 5]) ;1];
        c_transformed(:,2) = matrix \ [boundsSource([1 3 6]) ;1];
        c_transformed(:,3) = matrix \ [boundsSource([1 4 5]) ;1];
        c_transformed(:,4) = matrix \ [boundsSource([1 4 6]) ;1];
        c_transformed(:,5) = matrix \ [boundsSource([2 3 5]) ;1];
        c_transformed(:,6) = matrix \ [boundsSource([2 3 6]) ;1];
        c_transformed(:,7) = matrix \ [boundsSource([2 4 5]) ;1];
        c_transformed(:,8) = matrix \ [boundsSource([2 4 6]) ;1];
        
        
        % new bounds
        m=min(c_transformed,[],2);
        M=max(c_transformed,[],2);
        if (referenceImage==referenceImage)
            newSpacing = referenceImage.spacing;
        else
            newSpacing= inputImg.spacing;
        end
        newSize = ceil((M(1:3)-m(1:3))./newSpacing)+1;
        newOrigin = m(1:3);
        newOrigin = newOrigin - rem(newOrigin-inputImg.origin,newSpacing)-newSpacing;
        
        outputImage = ImageType(newSize,newOrigin,newSpacing,eye(3));
        
        
    else
         if (referenceImage==referenceImage)
            outputImage = ImageType(referenceImage);
         else
             outputImage = ImageType(inputImg);
         end
    end
    
    %% Divide in chunks to save memory
    
    
    NCHUNKS = ceil(outputImage.size/MAX_CHUNK_SIZE);
    
    
    chunked_size = ceil(outputImage.size./NCHUNKS)';
 
    [ix, iy, iz]= ndgrid(0:NCHUNKS(1)-1,0:NCHUNKS(2)-1,0:NCHUNKS(3)-1);
    intervals = [ix(:) iy(:) iz(:)];
    clear ix iy iz;
    for i=1:size(intervals,1)
         ranges([1 3 5]) = intervals(i,:).*chunked_size+1;
        ranges([2 4 6]) = min([(intervals(i,:)+[1 1 1]).*chunked_size ; outputImage.size']);
     
        ranges_size = ranges([2 4 6])-ranges([1 3 5])+[1 1 1];
        % generate all the indexes of the target image 
        [x, y, z] = ndgrid( ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6));
        
        % generate all the positions of the output image
        positionInTarget = outputImage.GetPosition([x(:) y(:) z(:)]');
     
         positionInSource = matrix * [positionInTarget ; ones( 1, size(positionInTarget,2) ) ];
         positionInSource(4,:)=[];

         valueInSource =  inputImg.GetValue(positionInSource,mode);
        % valueInSource(valueInSource~=valueInSource)=0;% remove NaN
        outputImage.data(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = reshape(valueInSource,ranges_size);
    end
     if (crop)
        outputImage = cropImage(outputImage, outputImage.GetBounds(0.99));% not 1 because this may produce issues with masks
     end

end
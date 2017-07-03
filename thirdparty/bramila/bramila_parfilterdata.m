function bramila_parfilterdata(cfg,s,mask,b,N,numfiles)
%PARFILTERDATA Subfunction of FILTERDATA
%   parfilterdata(cfg,s,mask,b,N,numfiles) contains the parallel iteration 
%   segment of the filterdata-function. Refer to filterdata.m for use of 
%   the filter.
%

% Display progress
if cfg.showprogress
    fprintf('Filtering file %1.0f / %1.0f\n',s,numfiles); 
end

runonce=1;

% Get input filename
if isfield(cfg,'infile')
    img = load_nii(cfg.infile{s});
    folders = strfind(cfg.infile{s},'/');
    filename = cfg.infile{s}(folders(end)+1:end);
else
    img = load_nii([cfg.inpath cfg.dirfile{s}]);
    filename = cfg.dirfile{s};
end

% Create output file
system(['touch ' cfg.outpath 'filtered_' filename]);

% Apply filter
for x=1:size(img.img,1);
    %disp(num2str(x)) % Display progress for testing
    for y=1:size(img.img,2);
        for z=1:size(img.img,3);
            if(mask.img(x,y,z)==0)
                continue
            end
            ts=double(squeeze(img.img(x,y,z,:)));
            me=mean(ts);
            ts=ts-me;
            temp=conv(ts,b);
            temp=flipud(conv(b,flipud(temp)));
            temp=temp+me;
            if(runonce==1) % Create blank image matrix
                L=length(temp((N+1):(end-N)));
                tempimg=zeros(...
                    size(img.img,1),...
                    size(img.img,2),...
                    size(img.img,3),L);
                runonce=0;
            end
            tempimg(x,y,z,:)=temp((N+1):(end-N));
        end
    end
end
img.img=tempimg;
img.hdr.dime.dim(5)=L;
img.hdr.dime.bitpix=64;
img.hdr.dime.datatype=64;

% Save output file
%fprintf('\nsave_nii(img,%s)\n',[cfg.outpath filename]);
save_nii(img,[cfg.outpath 'filtered_' filename]);

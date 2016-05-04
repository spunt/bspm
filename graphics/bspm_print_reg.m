function bspm_print_reg(in,prefix)

% ========================================================================%
if nargin<2, prefix = 'printreg'; end
if nargin<1, mfile_showhelp; return; end
if ischar(in), in = cellstr(in); end
allaxial = [];
allsagittal = [];
for i = 1:length(in)
    
    h = spm_vol(in{i});
    d = spm_read_vols(h);
    slices = round(size(d)/2);
    axial = d(:,:,slices(3));
    allaxial = [allaxial fliplr(axial)'];
    sagittal = squeeze(d(slices(1),:,:));
    allsagittal = [allsagittal fliplr(sagittal)'];
    
end
allaxial = rescale_im(allaxial);
allsagittal = rescale_im(allsagittal);
imwrite(allaxial,[prefix '_axial.jpg'],'jpg');
imwrite(allsagittal,[prefix '_sagittal.jpg'],'jpg');
function out = rescale_im(in)
out = in;
[r c] = size(in);
mn = min(in);
mx = max(in);
for i = 1:c
    out(:,i) = (in(:,i) - mn(i))/(mx(i)-mn(i));
end

    
        
        
    






 
 
 
 

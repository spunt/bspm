function out = matrixToRGB(input, cmap,limits)

ymin = 1;
ymax = size(cmap,1);

% if I do not specify otherwise...
xmin = limits(1);
xmax = limits(2);
s = size(input);

input_ = input(:);
indices  = min( [ max( [ round( (input_  -xmin)*(ymax-ymin)/(xmax-xmin) +ymin) ones(size(input_)) ],[],2) ymax*ones(size(input_))],[],2);

out_r = cmap(indices,1);
out_g = cmap(indices,2);
out_b = cmap(indices,3);

out = cat(numel(s)+1,reshape(out_r,s) , reshape(out_g,s), reshape(out_b,s));

end


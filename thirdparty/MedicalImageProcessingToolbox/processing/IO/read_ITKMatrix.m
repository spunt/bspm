function m =read_ITKMatrix(filename)

 fid = fopen(filename, 'r');
 C = textscan(fid, '%s', 4*4, 'HeaderLines', 1)';
 c = str2double(C{1});
 m = reshape(c,4,4)';   

end


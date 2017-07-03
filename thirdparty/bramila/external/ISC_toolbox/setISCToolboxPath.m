function setISCToolboxPath

p = path;
c = cd;
if isunix
    sl = '/';
else
    sl = '\';
end
c = [c sl];

path0 = cd;
path1 = [c 'ISCanalysis'];
path2 = [c 'visGUI'];
path3 = [c 'niftitools'];
path4 = [c 'templates'];

lastwarn(''); 
p = path(p,path1);
p = path(p,path2);
p = path(p,path3);
if exist(path4) == 7 % JT 8.11
  p = path(p,path4);
  flag = 1;
else 
  flag = 0;	
end	
if isempty(lastwarn)
    p = path(p,path0);
    disp('Added the following paths:')
    disp(path0)
    disp(path1)
    disp(path2)
    disp(path3)
    if flag  
      disp(path4)
    end  
end




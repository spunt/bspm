function m =read_CHeartMesh(filenameT,filenameX,varargin)
%

ncomp=3;
is_quad = false;
i=1;
while (i <= size(varargin,2))
    if (strcmp( varargin{i} , 'components'))
        ncomp= varargin{i+1};
        i = i+1;
    elseif (strcmp( varargin{i} , 'quadratic'))
        is_quad= true;
    end
    i = i+1;
end


fid = fopen(filenameX,'r');
M = fscanf(fid,'%i',[2 1]);
Nodes = fscanf(fid, '%f',[ncomp M(1)]);
fclose(fid);
Nodes = Nodes';

fid = fopen(filenameT,'r');
MT = fscanf(fid,'%i',[2 1]);
if (is_quad)
    % todo replace this for quadratic meshes
    T_ = fscanf(fid, '%i'); 
    T =reshape(T_,10,[]);
    T = T([1 2 3 4 5 7 6 8 9 10],:); % invert 6 and 7
%     T = [Tr(:,1:4) 
%          Tr(:, [1 5 7 8])
%          Tr(:, [2 6 5 9])
%          Tr(:, [2 7 6 10])
%          Tr(:, [8 9 10 4])]';
    
else
    %LETTER = fscanf(fid,'%s',[3 1]);
    T = fscanf(fid, '%i',[4 MT(1)]);
end
fclose(fid);
T = T';

m.Nodes = Nodes;
m.T= T;

end


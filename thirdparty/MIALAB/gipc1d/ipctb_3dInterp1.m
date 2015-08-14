function arOut = ipctb_3dInterp1(arIn, nExpand, nDim)
% nDim is the direction you want to expand the arIn-matrix nExpand times
% Cyrus 4/18/2008

nSizeX = size(arIn, 1);
nSizeY = size(arIn, 2);
nSizeZ = size(arIn, 3);
nExpX=1; nExpY=1; nExpZ=1;

switch nDim
    case 1, nMax1 = size(arIn, 2); nMax2 = size(arIn, 3);nExpX=nExpand*nSizeX-1;
    case 2, nMax1 = size(arIn, 1); nMax2 = size(arIn, 3);nExpY=nExpand*nSizeY-1;
    case 3, nMax1 = size(arIn, 1); nMax2 = size(arIn, 2);nExpZ=nExpand*nSizeZ-1;
    otherwise, error('Bug in ipctb_3dIntertp1');
end

arOut = zeros(nExpX, nExpY, nExpZ);
for i1 = 1:nMax1
    for i2 = 1:nMax2
        switch nDim
            case 1, arOut(:, i1, i2) = interp1((1:nSizeX)', arIn(:, i1, i2), linspace(1, nSizeX, nExpX)', 'linear');
            case 2, arOut(i1, :, i2) = interp1((1:nSizeY), arIn(i1, :, i2), linspace(1, nSizeY, nExpY), 'linear');
            case 3, arOut(i1, i2, :) = squeeze(interp1((1:nSizeZ), arIn(i1, i2, :), linspace(1, nSizeZ, nExpZ), 'linear'));
            otherwise, error('Bug in ipctb_3dIntertp1');
        end        
    end
end
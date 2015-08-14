function bRet = ipctb_sloverDisp(sType, sFileMap, sFileSave, rodTVals, hSlover)
% Exports nifti-files to slover and image file
global gIpc;

cmGray = [linspace(0,1, 64)' linspace(0,1, 64)' linspace(0,1, 64)'];

%spm_figure

bRet = false;
oSlover = slover;
oSlover.figure_struct.Units = 'pixels';
oSlover.slices = gIpc.ronSloverSlices;
oSlover.transform = eye(4,4);
oSlover.refreshf = 1;
oSlover.clf = 1;
oSlover.resurrectf = 1;
oSlover.userdata = 1;
oSlover.area.units='normalized';
oSlover.area.valign='top';
oSlover.area.halign='center';
oSlover.xslices = [];
oSlover.labels.colour = [1 1 1];
oSlover.labels.format = '%+3.0f';
oSlover.img(1).background = NaN;
oSlover.img(2).background = NaN;

switch(sType)
    case('blob')
        oSlover.cbar = 2;
        cmHot = [gIpc.sPathApp ipctb_backslash(gIpc.sMaskDir) 'cmHot.mat'];
        oSlover.img(1).type = 'truecolour';  
        oSlover.img(1).vol = ipctb_spm_vol(gIpc.sAnatTemp_dir);
        oSlover.img(1).cmap = cmGray;   
        oSlover.img(2).type = 'split';          
        oSlover.img(2).vol = ipctb_spm_vol(sFileMap);
        oSlover.img(2).cmap = cmHot;
        oSlover.img(2).outofrange = {0 64};        
        oSlover.img(2).range = rodTVals;        
    case('clust')
        oSlover.cbar = 2;        
        oSlover.img(1).type = 'truecolour';                
        oSlover.img(1).cmap = cmGray;           
        oSlover.img(1).vol = ipctb_spm_vol(gIpc.sAnatTemp_dir);
        oSlover.img(2).outofrange = {0 gIpc.nClusters};
        oSlover.img(2).range = [0.9 ; gIpc.nClusters+.1];        
        oSlover.img(2).cmap = gIpc.cmClust(2:gIpc.nClusters+1, :);
        oSlover.img(2).type = 'split';
        oSlover.img(2).hold = 0;        
        oSlover.img(2).vol = ipctb_spm_vol(sFileMap);              
    otherwise
        bRet = false;
        return
end

% Save image
oSlover.figure = hSlover;
oSlover = paint(oSlover);
% if isunix
%     eval(['print -f' num2str(oSlover.figure) ' -dpng ' '''' sFileSave ''''])
% else
set(oSlover.figure, 'PaperPositionMode', 'auto');
eval(['print -f' num2str(oSlover.figure) ' -dpng ' '''' sFileSave ''''])
%saveas(oSlover.figure, sFileSave);
% end
% Cyrus 5/7/2008 tried to save using the following
% but it often glitched in both pc and unix
% and saved the GIFT logo instead of brain
%     print_cmd='print -dpng -painters -noui';
%     print_fig(oSlover, sFileSave, print_cmd);
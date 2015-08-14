function PPPI_FSFAST_WRAPPER(subjectnos)

addpath('')
addpath('')
addpath('')

Study_Dir='';

Subjects={};

regionfile={};

regionname={};

models={'.lh' '.rh' '.mni305'};

for ii=subjectnos
    try
        clear config 
    catch
    end
    
    config.Subject=Subjects{ii};
    config.Study=Study_Dir;
    config.Subject_Subdirectory='bold';
    config.Model=''; % change for new model
    config.zip=0;
    
    try
        FSFAST2SPM(config)
    catch
        continue; %go to next subject
    end
    
    
    for jj=1:numel(regionfile)
        for mm=1:3
        %try
            Directory=[Study_Dir filesep Subjects{ii} filesep config.Subject_Subdirectory filesep 'SPMana' filesep config.Model models{mm} ];
            cd(Directory)
            load('/pathto/PPI_configfilename.mat');
            P.subject=Subjects{ii};
            P.directory=Directory;
            P.VOI=regionfile{jj};
            P.Region=regionname{jj};
            P.FSFAST=1;
            P.FSFAST_source=[Study_Dir filesep Subjects{ii} filesep config.Subject_Subdirectory filesep 'SPMana' filesep config.Model models{3} filesep 'SPM.mat']; % subcortical iamge for ROI extraction
            P.wb=1; % may need to be changed to 0 on local computers
            P.concatR=1;
            save([Subjects{ii} '_PPI_config.mat'],'P');
            PPPI([Subjects{ii} '_PPI_config.mat'])
        %catch
            disp(['Failed: ' Subjects{ii}])
        %end
        end
    end
    % You will want to zip the files here if P.zipfiles is not set to 1.
end


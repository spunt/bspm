function [cfg,group_cfg] = bramila_diagnostics(varargin)
%BRAMILA_ Summary of this function goes here
%   Detailed explanation goes here
% cfg.plot = 1; Generates diagnostic picture as in Power et al 2014

fprintf('\nData diagnostics\n');
cfg=varargin{1};
if nargin==1 % Single subject way
    
    tempcfg{1} = cfg;
    group_cfg = 0;
else
    tempcfg = cfg;
    group_cfg = varargin{2};
end

for subj = 1:length(tempcfg)
    fprintf('\nSubject %i (FileID ''%s'')\n',subj,tempcfg{subj}.fileID);
    % compute quality control measures
    tempcfg{subj}.dvars = bramila_dvars(tempcfg{subj});
    [tempcfg{subj}.fDisplacement,tempcfg{subj}.fDisplacement_rms]=bramila_framewiseDisplacement(tempcfg{subj});
    tempcfg{subj}.artrep_bad_volumes = bramila_artrep(tempcfg{subj});
    if(~isfield(tempcfg{subj},'plot') || (isfield(tempcfg{subj},'plot') && tempcfg{subj}.plot == 1) )
        disp('Plotting diagnostic images')
        % check that we have a mask and that size of the mask is compatible
        if(~isfield(tempcfg{subj},'mask'))
            warning('There is no mask for this subject, diagnostic plots will not be performed');
            continue;
        else
            if(all(size(tempcfg{subj}.mask)==[91 109 91]))
                % it only works for the 2mm MNI size
                atlasmasknii=load_nii('external/HO/mask.nii'); %1==WM, 2==CSF, 3==GM, 4==SUB, 5==CEREBELL
                atlasmask=atlasmasknii.img(:,:,:,1).*tempcfg{subj}.mask;
                vol=tempcfg{subj}.vol;
                T=size(vol,4);
                data=reshape(vol,[],T);
                data=data'; % Time in 1st dimension
                data=bramila_bold2perc(data);
                if(isfield(tempcfg{subj},'analysis_mask') && isfield(tempcfg{subj},'white_mask') && isfield(tempcfg{subj},'csf_mask'))
                    if(ischar(tempcfg{subj}.white_mask))
                        white_mask=load_nii(tempcfg{subj}.white_mask);
                    end
                    wmIDs=find(white_mask.img>tempcfg{subj}.white_mask_th);
                    if(ischar(tempcfg{subj}.csf_mask))
                        csf_mask=load_nii(tempcfg{subj}.csf_mask);
                    end
                    csfIDs=find(csf_mask.img>tempcfg{subj}.csf_mask_th);
                    if(ischar(tempcfg{subj}.analysis_mask))
                        analysis_mask_nii=load_nii(tempcfg{subj}.analysis_mask);
                    end
                    gmIDs=find(tempcfg{subj}.analysis_mask>0);
                else
                    disp('Masks were not computed, using a standard atlas...')
                    wmIDs=find(atlasmask==1);
                    csfIDs=find(atlasmask==2);
                    gmIDs=find(atlasmask>=3);
                end
                gsIDs=find(tempcfg{subj}.mask);
                GS=mean(data(:,gsIDs),2); % we could compute it from BOLD values and then percentages
                SD=std(data(:,gsIDs),0,2);
                DV=tempcfg{subj}.dvars;
                FD=tempcfg{subj}.fDisplacement;
                GM=mean(data(:,gmIDs),2);
                WM=mean(data(:,wmIDs),2);
                CSF=mean(data(:,csfIDs),2);
                
               
                
                figure
                % plot of FD
                subplot(4,1,1);
                bad5=100*length(find(FD>0.5))/T;
                bad2=100*length(find(FD>0.2))/T;
                plot(FD,'r');
                M=max([max(FD) 0.5]);
                axis([1 T 0 M]);
                text(5,M-M/100,['% Vols > 0.5: ' num2str(round(100*bad5)/100)],'VerticalAlignment','top');
                text(5,M*0.8-M/100,['% Vols > 0.2: ' num2str(round(100*bad2)/100)],'VerticalAlignment','top');
                ylabel('mm')
                dtemp=daspect;
                dtemp(1)=dtemp(1)/5;
                daspect(dtemp);
                ih=colorbar('EastOutside');
                set(ih,'Visible','off');
                legend('Framewise Displacement','Location','NorthEast')
                
                
                % plot of DV SD GS
                subplot(4,1,2);
                plot(DV,'b');
                hold on
                plot(SD,'g');
                plot(GS,'k');
                axis([1 T -5 5]);
                ylabel('BOLD');
                dtemp=daspect;
                dtemp(1)=dtemp(1)/5;
                daspect(dtemp);
                ih=colorbar('EastOutside');
                set(ih,'Visible','off');
                legend('DVARs','StdDev','Global Signal','Location','SouthEast','Orientation','Horizontal')
                
                % plot of some GM
                subplot(4,1,3);
                step=round(length(gmIDs)/800);
                ih=imagesc((data(:,gmIDs(1:step:end)))',[-5 5]);
                hold on
                text(5,5,'Grey matter','Color',[1 1 1],'VerticalAlignment','top');               
                colormap('gray');
                ylabel('GM');
                set(gca,'YTick',[]);
                ih=colorbar('EastOutside');
                ylabel(ih,'BOLD');
                dtemp=daspect;
                dtemp(1)=dtemp(1)/5;
                daspect(dtemp);
                
                
                subplot(4,1,4);
                step=round(length(wmIDs)/800);
                tempimg=[data(:,wmIDs(1:step:end))];
                WML=size(tempimg,2);
                step=round(length(csfIDs)/100);
                tempimg=[tempimg   data(:,csfIDs(1:step:end))];
                ih=imagesc(tempimg',[-5 5]);
                hold on
                plot([1 T],[WML WML],'Color',[1 1 1]);
                text(5,5,'White matter','Color',[1 1 1],'VerticalAlignment','top');
                text(5,WML+5,'CSF','Color',[1 1 1],'VerticalAlignment','top');
                ylabel('WM & CSF');
                set(gca,'YTick',[]);
                xlabel('Volume #');
                colormap('gray');
                ih=colorbar('EastOutside');
                ylabel(ih,'BOLD');
                dtemp=daspect;
                dtemp(1)=dtemp(1)/5;
                daspect(dtemp);
                saveas(gcf,[cfg.outpath '/bramila/diagnostics_time.eps'],'psc2');
                figure
                diagnosticMatrix=corr([ FD DV SD GS GM WM CSF]);
                dmlabels={'FD','DV','SD','GS','GM','WM','CSF'};
                imagesc(diagnosticMatrix,[-1 1])
                set(gca,'YTick',1:length(dmlabels))
                set(gca,'XTick',1:length(dmlabels))
                set(gca,'YtickLabel',dmlabels);
                set(gca,'XtickLabel',dmlabels);
                map=cbrewer('div','RdBu',15);
                map=flipud(map);
                colorbar
                colormap(map)
                
                saveas(gcf,[cfg.outpath '/bramila/diagnostics_corrMat.eps'],'psc2');
                save([cfg.outpath '/bramila/diagnostics.mat'],'gmIDs','wmIDs','csfIDs', 'gsIDs','FD', 'DV', 'SD', 'GS', 'GM', 'WM', 'CSF');
            end
        end
    end
end

if nargin==1
    cfg=tempcfg{1};
else
    cfg=tempcfg;
end
end


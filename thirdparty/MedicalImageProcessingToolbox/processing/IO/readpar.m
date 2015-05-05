function [par] = readpar(parfile)
%[par] = readpar(parfile);
%read scan parameters from exported PAR files (Philips research tools)
%
% output :
% - par   : scan parameters and image scaling factors - 
%           float_data = (int_data * RS + RI) / (RS*SS)
% 

par.problemreading=0; % will be set to 1 when an error occurs
parameter = textread (parfile,'%s',8, 'headerlines',7);
par.ResToolsVersion=parameter{8};

switch(par.ResToolsVersion)
    case 'V3'
        parameter = textread (parfile,'%s',5, 'delimiter', '.:','headerlines',11);
        par.name = parameter (4);
        parameter = textread (parfile,'%u', 16,'delimiter','.:Acquisitionr','headerlines',15);
        par.scno = parameter (16);
	parameter = textread (parfile,'%u', 28,'delimiter','.:Max.numberofcardiacphases','headerlines',18);
        par.phases = parameter (28);
	parameter = textread (parfile,'%u', 21,'delimiter','.:Max.numberofechoes','headerlines',19);
        par.echoes = parameter (21);
        parameter = textread (parfile,'%u', 31,'delimiter','.:Max.numberofslices/locations','headerlines',20);
        par.slice = parameter (31);
        parameter = textread (parfile,'%u', 23,'delimiter','.:Max.numberofdynamics','headerlines',21);
        par.dyn = parameter (23);
        parameter = textread (parfile,'%s', 25,'delimiter','.:Imagepixelsize[orbits]','headerlines',23);
        par.bit = parameter {25};
        parameter = textread (parfile,'%u', 24,'delimiter','.:Reconresolution(x,y)','headerlines',28);
        %parameter = textread (parfile,'%u', 24,'delimiter','.:Scanresolution(x,y)','headerlines',26);
        x = parameter (23);
        y = parameter (24);
        z = par.slice;
        par.dim = [x,y,z];
        parameter = textread (parfile,'%s', 20,'headerlines',90);
        par.sliceorient  = parameter (20);
        parameter = textread (parfile,'%f', 22,'delimiter','.:FOV(ap,fh,rl)[mm]','headerlines',31);

        if strcmp (par.sliceorient, '1')  ; % slice orientation: transversal
            fovx = parameter (20);
            fovy = parameter (22);
            par.sliceorient;
        end;

        if strcmp (par.sliceorient, '2')   ;% slice orientation: sagital
            fovx = parameter (21);
            fovy = parameter (20);
            par.sliceorient;
        end;

        if strcmp (par.sliceorient,'3')   ;% slice orientation: coronal
            fovx = parameter(22);
            fovy = parameter(21);
            par.sliceorient;
        end;

        parameter = textread (parfile,'%f', 21,'delimiter','.:Slicethickness[mm]','headerlines',32);
        par.slth = parameter (21);
        parameter = textread (parfile,'%f', 15,'delimiter','.:Slicegap[mm]','headerlines',33);
        par.gap= parameter (15);
        fovz = (par.gap + par.slth)*par.slice;
        par.fov = [fovx,fovy,fovz];
        parameter = textread (parfile,'%f', 39,'delimiter','.:Angulationmidslice(ap,fh,rl)[degr]','headerlines',35);
        par.angAP = parameter (37);
        par.angFH = parameter (38);
        par.angRL = parameter (39);
        parameter = textread (parfile,'%f', 36,'delimiter','.:OffCentremidslice(ap,fh,rl)[mm]','headerlines',36);
        par.offAP= parameter (34);
        par.offFH= parameter (35);
        par.offRL= parameter (36);
        parameter = textread (parfile,'%s',24, 'headerlines',88);
        voxx = str2num(parameter{23});
        voxy = str2num(parameter{24});
        voxz = par.slth + par.gap;
        par.vox=[voxx voxy voxz];
        parameternextline = textread (parfile,'%s',24, 'headerlines',90);
        if (parameternextline{1}-parameter{1})>0,
            par.slicessorted=1;
        else
            par.slicessorted=2;
        end

        clear parameter;
	
	% read scaling factors
	par.ri = []; 
	par.rs = []; 
	par.ss = [];
	par.sl = []; par.ec = []; par.dyn = []; par.ph = []; par.ty = [];
	par.echo_t = []; 
	par.tfe = [];
	fid = fopen(parfile);
	while (~feof(fid)),
	  str=fgetl(fid);
	  
	  [A,countA] = sscanf(str,' %d  %d   %d  %d %d %d   %d  %f %f %f  %d  %d  %f %f   %f  %f  %f   %f %d %d %d %d  %f  %f   %f   %f  %f     %f  %f');
	  if(~isempty(A) & (countA==29)), 
	    par.sl  = [par.sl; A(1)];  % slice number
	    par.ec  = [par.ec; A(2)];  % echo number
	    par.dyn = [par.dyn; A(3)];  % dynamic number
	    par.ph  = [par.ph; A(4)];  % phase number
	    par.ty  = [par.ty; A(5)];  % image type
	    par.echo_t = [par.echo_t; A(20)]; % echo time
	    par.ri = [par.ri; A(8)];
	    par.rs = [par.rs; A(9)];
	    par.ss = [par.ss; A(10)];
	    %par.tfe = [par.tfe; A(29)];
	  end
	end
	fclose(fid);
	
	% replace image_type : par.ty = 1 to max_number of types
	% (avoid having jumps in par.ty)
	% ---------------------------------------------
	k = 1;
	for l=0:max(par.ty),
	  idt = find(par.ty == l);
	  if ~isempty(idt),
	    par.ty(idt) = -k;
	    k = k+1;
	  end
	end
	par.ty = -par.ty;


    case 'V4'
        parameter = textread (parfile,'%s',5, 'delimiter', '.:','headerlines',11);
        par.name = parameter (3);
        parameter = textread (parfile,'%u', 16,'delimiter','.:Acquisitionr','headerlines',16);
        par.scno = parameter (16);
	parameter = textread (parfile,'%u', 28,'delimiter','.:Max.numberofcardiacphases','headerlines',19);
        par.phases = parameter (28);
	parameter = textread (parfile,'%u', 21,'delimiter','.:Max.numberofechoes','headerlines',20);
        par.echoes = parameter (21);
	parameter = textread (parfile,'%u', 31,'delimiter','.:Max.numberofslices/locations','headerlines',21);
        par.slice = parameter (31);
        parameter = textread (parfile,'%u', 23,'delimiter','.:Max.numberofdynamics','headerlines',22);
        par.dynamic = parameter (23);
	parameter = textread (parfile,'%u', 23,'delimiter','.:Scanresolution(x,y)','headerlines',27);
        par.scanres_x = parameter (22);
	par.scanres_y = parameter (23);
        % read first six rows of indices from per-slice lines in PAR file
        [slice_index(:,1) slice_index(:,2) slice_index(:,3) slice_index(:,4) slice_index(:,5) slice_index(:,6) ...
            slice_index(:,7) slice_index(:,8) slice_index(:,9) slice_index(:,10) slice_index(:,11)] = textread (parfile,'%n%n%n%n%n%n%n%n%n%n%n%*[^\n]','delimiter',' ','headerlines',90,'commentstyle','shell');
        par.slice_index=slice_index;
        parameter = textread (parfile,'%s',41, 'headerlines',91); % read first line of slice acquisition info
        if length(parameter)<15,
            disp(sprintf('Problem: No actual slices measured for volume %s.\nVolume might be corrupt.', parfile));
            par.problemreading=1;
        else
            par.sliceorient  = parameter{26};

            x = str2num(parameter{10});
            y = str2num(parameter{11});
            z = par.slice;
            par.dim = [x,y,z];

            par.bit = parameter{8};
            par.slth = str2num(parameter{23});
            par.gap = str2num(parameter{24});
            voxx = str2num(parameter{29});
            voxy = str2num(parameter{30});
            voxz = par.slth+par.gap;
            par.vox=[voxx voxy voxz];
            % to check whether slices are ordered (ie if all slices are
            % written sequentially, or all volumes)
            parameternextline = textread (parfile,'%s',24, 'headerlines',92);
            if (parameternextline{1}-parameter{1})>0,
                par.slicessorted=1;
            else
                par.slicessorted=2;
            end

            parameter = textread (parfile,'%f', 22,'delimiter','.:FOV(ap,fh,rl)[mm]','headerlines',30);

            if strcmp (par.sliceorient, '1')  ; % slice orientation: transversal
                fovx = parameter (20);
                fovy = parameter (22);
            end

            if strcmp (par.sliceorient, '2')   ;% slice orientation: sagital
                fovx = parameter (21);
                fovy = parameter (20);
            end

            if strcmp (par.sliceorient,'3')   ;% slice orientation: coronal
                fovx = parameter (22);
                fovy = parameter (21);
            end
            fovz = (par.gap + par.slth)*par.slice;
            par.fov = [fovx,fovy,fovz];


            parameter = textread (parfile,'%f', 36,'delimiter','.:Angulationmidslice(ap,fh,rl)[degr]','headerlines',32);
            par.angAP = parameter (34);
            par.angFH = parameter (35);
            par.angRL = parameter (36);
            parameter = textread (parfile,'%f', 36,'delimiter','.:OffCentremidslice(ap,fh,rl)[mm]','headerlines',33);




            par.offAP= parameter (34);
            par.offFH= parameter (35);
            par.offRL= parameter (36);

        end
	
	% read scaling factors
	par.ri = []; 
	par.rs = []; 
	par.ss = [];
	par.sl = []; par.ec = []; par.dyn = []; par.ph = []; par.ty = [];
	par.echo_t = []; 
	par.tfe = []; 
	par.reconres_x = []; par.reconres_y = []; % recon resolution (x y) 
	
	fid = fopen(parfile);
	while (~feof(fid)),
	  str=fgetl(fid);
	  %[A,countA] = sscanf(str,[' %d  %d   %d  %d %d %d   %d %d %d  %d %d  %f %f %f  %d  %d   %f   %f   %f   %f   %f %f  %f %f %f %d %d  %f %f   %f   %f %f  %f  %d  %f  %d  %d %d  %d   %d %d %d %d %f']);
	  [A,countA] = sscanf(str,[' %d  %d   %d  %d %d %d   %d %d %d  %d %d    %f  %f %f  %d  %d  %f   %f   %f   %f   %f   %f  %f %f  %d %d %d %d  %f %f   %f   %f %f  %f  %d  %f  %d %d %d %d   %f']);
	  if(~isempty(A) & (countA==41)),
	    par.sl  = [par.sl; A(1)];  % slice number
	    par.ec  = [par.ec; A(2)];  % echo number
	    par.dyn = [par.dyn; A(3)];  % dynamic number
	    par.ph  = [par.ph; A(4)];  % phase number
	    par.ty  = [par.ty; A(5)];  % image type
	    par.echo_t = [par.echo_t; A(31)]; % echo time
	    if isempty(par.reconres_x),
	      par.reconres_x = A(10);
	      par.reconres_y = A(11);
	    end
	    par.ri = [par.ri; A(12)];
	    par.rs = [par.rs; A(13)];
	    par.ss = [par.ss; A(14)];
	    par.tfe = [par.tfe; A(40)];
	  end
	end
	fclose(fid);
	
	% replace image_type : par.ty = 1 to max_number of types
	% (avoid having jumps in par.ty)
	% ---------------------------------------------
	k = 1;
	for l=0:max(par.ty),
	  idt = find(par.ty == l);
	  if ~isempty(idt),
	    par.ty(idt) = -k;
	    k = k+1;
	  end
	end
	par.ty = -par.ty;
    
    case 'V4.2'
        parameter = textread (parfile,'%s',5, 'delimiter', '.:','headerlines',11);
        par.name = parameter (3);
        %Ben: read in the protocol name
        parameter = textread (parfile,'%s',3,'delimiter','.:','headerlines',13);
        par.proname = parameter(3);
        %read in the kinetic time;
        tmp = textread (parfile,'%s',6,'headerlines',13);
        if strcmp(tmp(5),'Kinetic') || strcmp(tmp(5),'TFE-Kinetic') || strcmp(tmp(5),'kinetic')
            ktime = char(tmp(6));
            offc = find(ktime == 'm');
            par.ktime = str2num(ktime(1:offc-1));
        else
            par.ktime = [];
        end
        %Ben: read in scan date
        parameter = textread (parfile,'%s',22,'delimiter','Examinationdate/time','headerlines',14) ;
        par.scandate = parameter(21);
        %ben: read in scan time
        par.scantime = parameter(22);        
        parameter = textread (parfile,'%u', 16,'delimiter','.:Acquisitionr','headerlines',16);
        par.scno = parameter (16);
        %Ben: read in the scan duration.
        parameter = textread (parfile,'%f',20,'delimiter','.:ScanDuration[sec]','headerlines',18);
        par.scandur = parameter(20);
        parameter = textread (parfile,'%s',5, 'delimiter', '.:','headerlines',13);
        par.prot = parameter (3);
        parameter = textread (parfile,'%s',5, 'delimiter', ':/','headerlines',14);
        par.date = parameter(3);
        parameter = textread (parfile,'%u', 16,'delimiter','.:Acquisitionr','headerlines',16);
        par.scno = parameter (16);
	parameter = textread (parfile,'%u', 28,'delimiter','.:Max.numberofcardiacphases','headerlines',19);
        par.phases = parameter (28);
	parameter = textread (parfile,'%u', 21,'delimiter','.:Max.numberofechoes','headerlines',20);
        par.echoes = parameter (21);
	parameter = textread (parfile,'%u', 31,'delimiter','.:Max.numberofslices/locations','headerlines',21);
        par.slice = parameter (31);
        parameter = textread (parfile,'%u', 23,'delimiter','.:Max.numberofdynamics','headerlines',22);
        par.dynamic = parameter (23);
	parameter = textread (parfile,'%u', 23,'delimiter','.:Scanresolution(x,y)','headerlines',27);
        par.scanres_x = parameter (22);
	par.scanres_y = parameter (23);
        % read first six rows of indices from per-slice lines in PAR file
        [slice_index(:,1) slice_index(:,2) slice_index(:,3) slice_index(:,4) slice_index(:,5) slice_index(:,6) ...
            slice_index(:,7) slice_index(:,8) slice_index(:,9) slice_index(:,10) slice_index(:,11)] = textread (parfile,'%n%n%n%n%n%n%n%n%n%n%n%*[^\n]','delimiter',' ','headerlines',97,'commentstyle','shell');
        par.slice_index=slice_index;
        parameter = textread (parfile,'%s',41, 'headerlines',99); % read first line of slice acquisition info
        if length(parameter)<15,
            disp(sprintf('Problem: No actual slices measured for volume %s.\nVolume might be corrupt.', parfile));
            par.problemreading=1;
        else
            par.sliceorient  = parameter{26};

            x = str2num(parameter{10});
            y = str2num(parameter{11});
            z = par.slice;
            par.dim = [x,y,z];

            par.bit = parameter{8};
            par.slth = str2num(parameter{23});
            par.gap = str2num(parameter{24});
            voxx = str2num(parameter{29});
            voxy = str2num(parameter{30});
            voxz = par.slth+par.gap;
            par.vox=[voxx voxy voxz];
            % to check whether slices are ordered (ie if all slices are
            % written sequentially, or all volumes)
            parameternextline = textread (parfile,'%s',24, 'headerlines',92);
            if (parameternextline{1}-parameter{1})>0,
                par.slicessorted=1;
            else
                par.slicessorted=2;
            end

            parameter = textread (parfile,'%f', 22,'delimiter','.:FOV(ap,fh,rl)[mm]','headerlines',30);

            if strcmp (par.sliceorient, '1')  ; % slice orientation: transversal
                fovx = parameter (20);
                fovy = parameter (22);
            end

            if strcmp (par.sliceorient, '2')   ;% slice orientation: sagital
                fovx = parameter (21);
                fovy = parameter (20);
            end

            if strcmp (par.sliceorient,'3')   ;% slice orientation: coronal
                fovx = parameter (22);
                fovy = parameter (21);
            end
            fovz = (par.gap + par.slth)*par.slice;
            par.fov = [fovx,fovy,fovz];


            parameter = textread (parfile,'%f', 39,'delimiter','.:Angulationmidslice(ap,fh,rl)[degr]','headerlines',32);
            par.angAP = parameter (37);
            par.angFH = parameter (38);
            par.angRL = parameter (39);
            parameter = textread (parfile,'%f', 36,'delimiter','.:OffCentremidslice(ap,fh,rl)[mm]','headerlines',33);
          
            par.offAP= parameter (34);
            par.offFH= parameter (35);
            par.offRL= parameter (36);
            
            parameter = textread(parfile,'%s',33,'delimiter',' :Phaseencodingvelocity[cm/sec]','headerlines',36);
            
            par.phaseEncX = str2double( parameter(31));
            par.phaseEncY = str2double(parameter(32));
            par.phaseEncZ = str2double(parameter(33));
            

        end
	
	% read scaling factors
	par.ri = []; 
	par.rs = []; 
	par.ss = [];
	par.sl = []; par.ec = []; par.dyn = []; par.ph = []; par.ty = []; par.sseq = [];
	par.echo_t = []; 
	par.tfe = []; 
	par.reconres_x = []; par.reconres_y = []; % recon resolution (x y) 
	par.freq = []; % cardiac frequency
    par.delay= []; % cardiac frequency
    par.spacing = [];
    par.slicetime = [];
	fid = fopen(parfile);
	while (~feof(fid)),
	  str=fgetl(fid);
	  [A,countA] = sscanf(str,['  %d   %d    %d  %d %d %d     %d  %d    %d  %d  %d   %f   %f %f   %d  %d   %f   %f   %f   %f   %f   %f  %f %f %d %d %d %d  %f  %f  %f    %f  %f    %f   %d   %f    %d  %d  %d   %d   %f %d    %d    %d    %d    %f     %f     %f']);
	  if(~isempty(A) & (countA==49)),
	    par.sl  = [par.sl; A(1)];  % slice number
	    par.ec  = [par.ec; A(2)];  % echo number
	    par.dyn = [par.dyn; A(3)];  % dynamic number
	    par.ph  = [par.ph; A(4)];  % phase number
	    par.ty  = [par.ty; A(5)];  % image type
        par.sseq  = [par.sseq; A(6)];  % image type
	    par.echo_t = [par.echo_t; A(31)]; % echo time
	    if isempty(par.reconres_x),
	      par.reconres_x = A(10);
	      par.reconres_y = A(11);
	    end
	    par.ri = [par.ri; A(12)];
	    par.rs = [par.rs; A(13)];
	    par.ss = [par.ss; A(14)];
	    par.tfe = [par.tfe; A(40)];
        par.freq = [par.freq  ;A(37)];
        par.spacing = [par.spacing ; A([29 30 23])];
        par.slicetime= [par.slicetime; A(33)];
        par.delay = [par.delay; A(32)];
	  end
	end
	fclose(fid);
	
	% replace image_type : par.ty = 1 to max_number of types
	% (avoid having jumps in par.ty)
	% ---------------------------------------------
    %  	k = 1;
    %  	for l=0:max(par.ty),
    %  	  idt = find(par.ty == l);
    %  	  if ~isempty(idt),
    %  	    par.ty(idt) = -k;
    %  	    k = k+1;
    %  	  end
    %  	end
    %  	par.ty = -par.ty;
    
    differentPairs=[];
    for i=1:numel(par.ty)
        currentPair = [par.ty(i) par.sseq(i)];
        [r c] =find(differentPairs==par.ty(i));
        [r2 c2]=find(differentPairs==par.sseq(i));
        
        
        if (numel(r)==0 || numel(r2)==0)
            differentPairs=[differentPairs; currentPair];
        end
        
    end

   
   % par.ty= size(differentPairs,1);

 
    case 'V4.1' 
        parameter = textread (parfile,'%s',5, 'delimiter', '.:','headerlines',11);
        par.name = parameter (3);
        %Ben: read in the protocol name
        parameter = textread (parfile,'%s',3,'delimiter','.:','headerlines',13);
        par.proname = parameter(3);
        %read in the kinetic time;
        tmp = textread (parfile,'%s',6,'headerlines',13);
        if strcmp(tmp(5),'Kinetic') || strcmp(tmp(5),'TFE-Kinetic') || strcmp(tmp(5),'kinetic')
            ktime = char(tmp(6));
            offc = find(ktime == 'm');
            par.ktime = str2num(ktime(1:offc-1));
        else
            par.ktime = [];
        end
        %Ben: read in scan date
        parameter = textread (parfile,'%s',22,'delimiter','Examinationdate/time','headerlines',14) ;
        par.scandate = parameter(21);
        %ben: read in scan time
        par.scantime = parameter(22);        
        parameter = textread (parfile,'%u', 16,'delimiter','.:Acquisitionr','headerlines',16);
        par.scno = parameter (16);
        %Ben: read in the scan duration.
        parameter = textread (parfile,'%f',20,'delimiter','.:ScanDuration[sec]','headerlines',18);
        par.scandur = parameter(20);
	parameter = textread (parfile,'%u', 28,'delimiter','.:Max.numberofcardiacphases','headerlines',19);
        par.phases = parameter (28);
	parameter = textread (parfile,'%u', 21,'delimiter','.:Max.numberofechoes','headerlines',20);
        par.echoes = parameter (21);
	parameter = textread (parfile,'%u', 31,'delimiter','.:Max.numberofslices/locations','headerlines',21);
        par.slice = parameter (31);
        parameter = textread (parfile,'%u', 23,'delimiter','.:Max.numberofdynamics','headerlines',22);
        par.dynamic = parameter (23);
	parameter = textread (parfile,'%u', 23,'delimiter','.:Scanresolution(x,y)','headerlines',27);
        par.scanres_x = parameter (22);
	par.scanres_y = parameter (23);
        % read first six rows of indices from per-slice lines in PAR file
        [slice_index(:,1) slice_index(:,2) slice_index(:,3) slice_index(:,4) slice_index(:,5) slice_index(:,6) ...
            slice_index(:,7) slice_index(:,8) slice_index(:,9) slice_index(:,10) slice_index(:,11)] = textread (parfile,'%n%n%n%n%n%n%n%n%n%n%n%*[^\n]','delimiter',' ','headerlines',97,'commentstyle','shell');
        par.slice_index=slice_index;
        parameter = textread (parfile,'%s',41, 'headerlines',99); % read first line of slice acquisition info
        if length(parameter)<15,
            disp(sprintf('Problem: No actual slices measured for volume %s.\nVolume might be corrupt.', parfile));
            par.problemreading=1;
        else
            par.sliceorient  = parameter{26};

            x = str2num(parameter{10});
            y = str2num(parameter{11});
            z = par.slice;
            par.dim = [x,y,z];

            par.bit = parameter{8};
            par.slth = str2num(parameter{23});
            par.gap = str2num(parameter{24});
            voxx = str2num(parameter{29});
            voxy = str2num(parameter{30});
            voxz = par.slth+par.gap;
            par.vox=[voxx voxy voxz];
            % to check whether slices are ordered (ie if all slices are
            % written sequentially, or all volumes)
            parameternextline = textread (parfile,'%s',24, 'headerlines',92);
            if (parameternextline{1}-parameter{1})>0,
                par.slicessorted=1;
            else
                par.slicessorted=2;
            end

            parameter = textread (parfile,'%f', 22,'delimiter','.:FOV(ap,fh,rl)[mm]','headerlines',30);
            if strcmp (par.sliceorient, '1')  ; % slice orientation: transversal
                fovx = parameter (20);
                fovy = parameter (22);
            end

            if strcmp (par.sliceorient, '2')   ;% slice orientation: sagital
                fovx = parameter (21);
                fovy = parameter (20);
            end

            if strcmp (par.sliceorient,'3')   ;% slice orientation: coronal
                fovx = parameter (22);
                fovy = parameter (21);
            end
            fovz = (par.gap + par.slth)*par.slice;
            par.fov = [fovx,fovy,fovz];


            parameter = textread (parfile,'%f', 36,'delimiter','.:Angulationmidslice(ap,fh,rl)[degr]','headerlines',32);
            par.angAP = parameter (34);
            par.angFH = parameter (35);
            par.angRL = parameter (36);
            parameter = textread (parfile,'%f', 36,'delimiter','.:OffCentremidslice(ap,fh,rl)[mm]','headerlines',33);




            par.offAP= parameter (34);
            par.offFH= parameter (35);
            par.offRL= parameter (36);
            
             parameter = textread(parfile,'%s',33,'delimiter',' :Phaseencodingvelocity[cm/sec]','headerlines',36);
            
            par.phaseEncX = str2double( parameter(31));
            par.phaseEncY = str2double(parameter(32));
            par.phaseEncZ = str2double(parameter(33));

        end
	
	% read scaling factors
	par.ri = []; 
	par.rs = []; 
	par.ss = [];
	par.sl = []; par.ec = []; par.dyn = []; par.ph = []; par.ty = [];
	par.echo_t = []; 
	par.tfe = []; 
	par.reconres_x = []; par.reconres_y = []; % recon resolution (x y) 
	
	fid = fopen(parfile);
    while (~feof(fid)),
	  str=fgetl(fid);
	 [A,countA] = sscanf(str,['  %d   %d    %d  %d %d %d     %d  %d    %f  %d  %d   %f   %f %f   %d  %d   %f   %f   %f   %f   %f   %f  %f %f %d %d %d %d  %f  %f  %f    %f  %f    %f   %f   %f    %d  %d  %d   %d   %f %d    %d    %d    %d    %f     %f     %f']);
      
      if(~isempty(A) && (countA==49)),
	    par.sl  = [par.sl; A(1)];  % slice number
	    par.ec  = [par.ec; A(2)];  % echo number
	    par.dyn = [par.dyn; A(3)];  % dynamic number
	    par.ph  = [par.ph; A(4)];  % phase number
	    par.ty  = [par.ty; A(5)];  % image type
	    par.echo_t = [par.echo_t; A(31)]; % echo time
	    if isempty(par.reconres_x),
	      par.reconres_x = A(10);
	      par.reconres_y = A(11);
	    end
	    par.ri = [par.ri; A(12)];
	    par.rs = [par.rs; A(13)];
	    par.ss = [par.ss; A(14)];
	    par.tfe = [par.tfe; A(40)];
      end
    end
    
	fclose(fid);
	
	% replace image_type : par.ty = 1 to max_number of types
	% (avoid having jumps in par.ty)
	% ---------------------------------------------
	k = 1;
	for l=0:max(par.ty),
	  idt = find(par.ty == l);
	  if ~isempty(idt),
	    par.ty(idt) = -k;
	    k = k+1;
	  end
	end
	par.ty = -par.ty;
end  

end





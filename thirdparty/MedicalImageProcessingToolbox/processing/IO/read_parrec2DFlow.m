function [outphase, outanatomy, patientdata] = read_parrec2DFlow(filename)


par_survey = parread(filename );

rec_survey = readrec([filename(1:end-3) 'REC' ]);

% Split the 5D 2Dflow data
data_mag = rec_survey(:,:,:,1,1);
data_oth = rec_survey(:,:,:,2,1); % ?
data_phs = rec_survey(:,:,:,2,2);

% calculate velocity values
max_ph = max(data_phs(:));
p = nextpow2(max_ph);
venc = max(par_survey.PhaseEncodingVelocity);
data_phs = data_phs - (2^p)/2;
data_phs = data_phs./((2^p)/2);
data_phs = data_phs.*(venc); % (division by 100 because cm/s?)


% Orientantion

AFRtoLPS = [0 0 1; 1 0 0; 0 1 0];

magicmatrix = [-1 0 0; 0 0 1; 0 -1 0];
TRA = [0 1 0; -1 0 0; 0 0 -1];
SAG = [-1 0 0; 0 0 1; 0 -1 0];
COR = [0 1 0; 0 0 1; -1 0 0];

switch (par_survey.ImageInformation.SliceOrientation(1))
    case 1
        Torientation = TRA;
    case 2
        Torientation = SAG;
    case 3
        Torientation = COR;
    otherwise
        Torientation = eye(3);
end

ap = par_survey.AngulationMidslice(1)* pi / 180.0;
fh = par_survey.AngulationMidslice(2)* pi / 180.0;
rl = par_survey.AngulationMidslice(3)* pi / 180.0;

Tap = [1 0 0; 0 cos(ap) -sin(ap); 0 sin(ap) cos(ap)];
Tfh = [cos(fh) 0 sin(fh); 0 1 0; -sin(fh) 0 cos(fh)];
Trl = [cos(rl) -sin(rl) 0; sin(rl) cos(rl) 0; 0 0 1];

orientation_matrix = AFRtoLPS * Trl * Tap * Tfh * magicmatrix' * Torientation'; % c++ version

s = size(data_phs)';
s_ = [s(1:2) ;1 ;s(3)];
% ORIGIN ----------------------------------------------------
tmax = 60/par_survey.ImageInformation.CardiacFrequency(1);
%extent = AFRtoLPS * par_survey.FOV';
extent = abs(Torientation * magicmatrix * par_survey.FOV'); % I am not sure exactly why I have to use the abs
spacing = [extent ; tmax]./s_;

midoffset =    spacing(1:3) .* (s(1:3)- ones(3,1)) ./2.0;
midoffset = orientation_matrix * midoffset;

origin = par_survey.OffCentreMidslice';
origin = AFRtoLPS * origin;

origin = origin - midoffset;

M = eye(4);
M(1:3,1:3)=orientation_matrix;


outphase = ImageType(s_,[origin ; 0],spacing,M);
outphase.data = reshape(data_phs, s_');
outanatomy = ImageType(s_,[origin ; 0],spacing,M);
outanatomy.data = reshape(data_mag, s_');


patientdata.type = '3DCDI';
%params.type = 'other';
patientdata.heart_rate = par_survey.ImageInformation.CardiacFrequency(1);
patientdata.frame_rate=s_(end)/(60/par_survey.ImageInformation.CardiacFrequency(1));
patientdata.trigger_delay = 0;
patientdata.venc = venc;

end

function values = parread(filename)

%PARREAD reads Philips .par file
%   PARREAD(FILENAME) reads the .par file and returns a struct
%   whose fields contain all parameters. PARREAD can handle all
%   versions of .par files
%
%   The field ImageInformation in the output struct is a struct
%   itself and contains the parameters which correspond to the
%   lower part of the .par file.
%
%   The fields of the output struct can be accessed by a .
%   Example:
%                   values = parread('My_Parfile')
%   reads the .par file and returns the output struct
%
%                   values.RepetitionTime
%   returns the repetition time
%
%                   values.ImageInformation.SliceNumber
%   returns the slice order in the .rec file

fid = fopen(filename);
if fid == -1
    filename = [filename(1:end-3),'par'];
    fid = fopen(filename);
    if fid == -1
        filename = [filename(1:end-3),'PAR'];
        fid = fopen(filename);
        if fid == -1
            error(['The .par file ',filename,' does not exist'])
            return
        end
    end
end
values = read_upper_part(fid);
values = read_lower_part(fid, values);
fclose(fid);
end
%--------- function for reading the upper part of the .par file ----------
function values = read_upper_part(fid)
%-------------------------------------------------------------------------
parameters = read_parameters(fid);
struct_names = format_parameters(parameters);
no_parameters = size(parameters,1);
values = struct;
h = waitbar(0,'Read header I');
for i = 1:no_parameters
    waitbar(i/no_parameters);
    values_temp = [];
    fseek(fid,0,-1);
    while ~feof(fid)
        
        s = fgetl(fid);
        s_index = findstr(parameters{i},s);
        if ~isempty(s_index)
            s_index = findstr(':', s);
            %             s_index = s_index(1)+4;
            s_index = s_index(1)+1;
            
            if isempty(str2num(s(s_index:length(s))))
                values_temp = s(s_index:length(s));
            else
                values_temp = str2num(s(s_index:length(s)));
            end
            values = setfield(values, struct_names{i},values_temp);
            s_index=s_index+1;
        end
    end
end
close(h);
end
%--------- function for reading the lower part of the .par file ----------
function values = read_lower_part(fid, values)
%-------------------------------------------------------------------------
[image_information_legend,image_information_length] = read_image_information_parameters(fid);
if ~isempty(image_information_legend)
    image_information_legend = format_parameters(image_information_legend);
    fseek(fid,0,-1);
    loop =1;
    while ~feof(fid)
        s = fgetl(fid);
        s_index = strmatch('# === IMAGE INFORMATION ==',s);
        if ~isempty(s_index)
            fgetl(fid);
            fgetl(fid);
            while ~feof(fid)
                s = fgetl(fid);
                if length(str2num(s)) ~= 0
                    h(loop,:) = str2num(s);
                end
                loop = loop+1;
            end
        end
    end
    if ~isempty(s_index)
        ImageInformation = struct;
        k = 1;
        for i = 1:length(image_information_legend)
            ImageInformation = setfield(ImageInformation,image_information_legend{i},h(:,k:k+image_information_length(i)-1));
            k = k+image_information_length(i);
        end
        values = setfield(values, 'ImageInformation',ImageInformation);
    end
end
end

%-- Check which parameters are stored in the upper part of the .par file --
function parameters = read_parameters(fid);
%--------------------------------------------------------------------------
fseek(fid,0,-1);
loop = 1;
while ~feof(fid)
    s = fgetl(fid);
    s_index = strmatch('# === GENERAL INFORMATION',s);
    if ~isempty(s_index)
        fgetl(fid);
        while ~feof(fid)
            s = fgetl(fid);
            if isempty(s) | ~strcmp(s(1),'.')
                break
            end
            ind = findstr(s,':');
            s = s(2:ind-1);
            parameters{loop} = strtrim(s);
            loop = loop +1;
        end
    end
end
parameters = parameters';
end
%-- Format parameter names such that they can be given as struct fields --
function struct_names = format_parameters(parameters)
%--------------------------------------------------------------------------
h = waitbar(0,'Read header II');
for i = 1:length(parameters)
    waitbar(i/length(parameters));
    s = parameters{i};
    ind1 = strfind(s,'(');
    ind2 = strfind(s,'[');
    ind3 = strfind(s,'<');
    ind = [ind1,ind2,ind3];
    if ~isempty(ind)
        ind = min(ind);
        s = s(1:ind-1);
    end
    s = strrep(s,'.',' ');
    s = strrep(s,'/',' ');
    s = strrep(s,'_',' ');
    s = strrep(s,'-',' ');
    s = deblank(s);
    ind = strfind(s,' ');
    s(ind+1) = upper(s(ind+1));
    s(1) = upper(s(1));
    s(ind)=[];
    struct_names{i} = s;
end
close(h);
struct_names = struct_names';
end
%-- Check which parameters are stored in the lower part of the .par file --
function [image_information_legend,image_information_length] = read_image_information_parameters(fid);
%--------------------------------------------------------------------------
fseek(fid,0,-1);
loop = 1;
image_information_legend = [];
image_information_length = [];
while ~feof(fid)
    s = fgetl(fid);
    s_index = strmatch('# === IMAGE INFORMATION DEFINITION',s);
    if ~isempty(s_index)
        fgetl(fid);
        fgetl(fid);
        while ~feof(fid)
            s = fgetl(fid);
            if length(s) < 5
                break
            end
            ind = max(findstr(s,'('));
            l = str2double(s(ind+1));
            s = s(2:ind-1);
            image_information_legend{loop} = strtrim(s);
            if isnan(l) | imag(l)~=0
                image_information_length(loop) = 1;
            else
                image_information_length(loop) = l;
            end
            loop = loop +1;
        end
    end
end
if ~isempty(image_information_legend)
    image_information_legend = image_information_legend';
end
end

function data = readrec(fn, read_params, v)

%------------------------------------------------------
% fReadrec: 	reads a rec-file
%
% Input:		fn				name of the rec-file
%               read_params     optional input. Specifies the images to be
%                               read. read_params is a struct created
%                               by the function create_read_param_struct
%               v               Parameter struct created by parread
%
% Output:	out_pics		dataset as a 3D-matrix
%------------------------------------------------------

parfile = [fn(1:end-3),'PAR'];

if nargin == 1
    [read_params, v, DataFormat] = ReadParameterFile(fn);
elseif nargin == 2
    v = parread(parfile);
end


if isfield(v,'ReconResolution')
    size1 = v.ReconResolution(1);
    size2 = v.ReconResolution(2);
else
    size1 = v.ImageInformation.ReconResolution(1);
    size2 = v.ImageInformation.ReconResolution(2);
end

fid=fopen(fn,'r','l');
% l means little ended

data = zeros(size1,size2,length(read_params.kz),length(read_params.echo),length(read_params.dyn),length(read_params.card),length(read_params.typ), length(read_params.mix),'single');
h = waitbar(0,'Read image data');
for sl = 1:length(read_params.kz)
    waitbar(sl/length(read_params.kz));
    for ec = 1:length(read_params.echo)
        for dy = 1:length(read_params.dyn)
            for ph = 1:length(read_params.card)
                for ty = 1:length(read_params.typ)
                    for mi = 1:length(read_params.mix)
                        ind{1} = find(v.ImageInformation.SliceNumber == read_params.kz(sl));
                        ind{2} = find(v.ImageInformation.EchoNumber == read_params.echo(ec));
                        ind{3} = find(v.ImageInformation.DynamicScanNumber == read_params.dyn(dy));
                        ind{4} = find(v.ImageInformation.CardiacPhaseNumber == read_params.card(ph));
                        ind{5} = find(v.ImageInformation.ImageTypeMr == read_params.typ(ty));
                        ind{6} = find(v.ImageInformation.ScanningSequence == read_params.mix(mi));
                        im_ind = ind{1};
                        for i = 2:6
                            im_ind = intersect(im_ind,ind{i});
                        end
                        offset = v.ImageInformation.IndexInRECFile(im_ind)*size1*size2*2;
                        if isempty(offset)
                            continue
                            %                             error('Some parameters specified are not valid')
                        end
                        fseek(fid,offset,-1);
                        im = fread(fid,size1*size2,'uint16');
                        data(:,:,sl,ec,dy,ph,ty,mi) = reshape(im,size1,size2)';
                    end
                end
            end
        end
    end
end
close(h);
data = squeeze(data);

end

function [Parameter2read, Parameter, DataFormat, Filenames] = ReadParameterFile(file)

dotind = findstr(file,'.');
if ~isempty(dotind)
    dotind = dotind(end);
    ending = lower(file(dotind(end)+1:end));
else
    ending = 'raw';
end

switch ending
    case {'rec', 'par'}
        parfile = [file(1:dotind),'par'];
        Parameter = parread(parfile);
        DataFormat = 'Rec';
        Filenames.DataFile = [file(1:dotind),'rec'];
        Filenames.ParameterFile = parfile;
        
        Parameter2read.kz = unique(Parameter.ImageInformation.SliceNumber);
        Parameter2read.echo = unique(Parameter.ImageInformation.EchoNumber);
        Parameter2read.dyn = unique(Parameter.ImageInformation.DynamicScanNumber);
        Parameter2read.card =  unique(Parameter.ImageInformation.CardiacPhaseNumber);
        Parameter2read.typ = unique(Parameter.ImageInformation.ImageTypeMr);
        Parameter2read.mix = unique(Parameter.ImageInformation.ScanningSequence);
        Parameter2read.aver = [];
        Parameter2read.rtop = [];
        Parameter2read.ky = [];
        Parameter2read.loca = [];
        Parameter2read.chan = [];
        Parameter2read.extr1 = [];
        Parameter2read.extr2 = [];
    case 'cpx'
        Parameter = read_cpx_header(file,'no');
        DataFormat = 'Cpx';
        Filenames.DataFile = file;
        Filenames.ParameterFile = file;
        
        Parameter2read.loca = unique(Parameter(:,1))+1;
        Parameter2read.kz = unique(Parameter(:,2))+1;
        Parameter2read.chan = unique(Parameter(:,3))+1;
        Parameter2read.card = unique(Parameter(:,4))+1;
        Parameter2read.echo = unique(Parameter(:,5))+1;
        Parameter2read.dyn = unique(Parameter(:,6))+1;
        Parameter2read.extr1 = unique(Parameter(:,7))+1;
        Parameter2read.extr2 = unique(Parameter(:,18))+1;
        Parameter2read.aver = [];
        Parameter2read.rtop = [];
        Parameter2read.typ = [];
        Parameter2read.mix = [];
        Parameter2read.ky = [];
    case {'data', 'list'}
        t = 'TEHROA';
        DataFormat = 'ExportedRaw';
        listfile = [file(1:dotind),'list'];
        Parameter = listread(listfile);
        
        Filenames.DataFile = [file(1:dotind),'data'];
        Filenames.ParameterFile = listfile;
        
        typ = unique(Parameter.Index.typ(:,2));
        for i = 1:length(typ)
            numtyp(i) = findstr(typ(i),t);
        end
        Parameter2read.typ = sort(numtyp);
        Parameter2read.mix = unique(Parameter.Index.mix);
        Parameter2read.dyn = unique(Parameter.Index.dyn);
        Parameter2read.card = unique(Parameter.Index.card);
        Parameter2read.echo = unique(Parameter.Index.echo);
        Parameter2read.loca = unique(Parameter.Index.loca);
        Parameter2read.chan = unique(Parameter.Index.chan);
        Parameter2read.extr1 = unique(Parameter.Index.extr1);
        Parameter2read.extr2 = unique(Parameter.Index.extr2);
        Parameter2read.ky = unique(Parameter.Index.ky);
        Parameter2read.kz = unique(Parameter.Index.kz);
        Parameter2read.aver = unique(Parameter.Index.aver);
        Parameter2read.rtop = unique(Parameter.Index.rtop);
    case {'raw', 'lab'}
        t = 'TEHROA';
        DataFormat = 'Raw';
        
        if isempty( dotind )
            slash_ind = strfind( file, filesep )+1;
            if isempty( slash_ind )
                slash_ind = 1;
                directory  = '';
                files = dir ;
            else
                directory = file( 1:slash_ind(end)-1 );
                files = dir( directory ) ;
            end
            one_file = file( slash_ind(end):end );
            
            ind_one_file = -1;
            ind_other_file = -1;
            for i = 1:size(files,1)
                if strcmp( files(i).name , one_file )
                    ind_one_file = i;
                end
                if strfind( files(i).name, one_file( 1:17)  ) & ...
                        isempty( strfind( files(i).name, '.' )) & ...
                        i ~= ind_one_file
                    ind_other_file = i;
                end                
            end
            if ind_one_file > 0 & ind_other_file > 0
                if files(ind_one_file).bytes > files(ind_other_file).bytes
                    Filenames.DataFile = [directory, files(ind_one_file).name];
                    Filenames.ParameterFile = [directory, files(ind_other_file).name];
                else
                    Filenames.DataFile = [directory,files(ind_other_file).name];
                    Filenames.ParameterFile = [directory,files(ind_one_file).name];
                end
            else
                error('data/parameter-file pair not found');
            end
            dotind = length(Filenames.ParameterFile)+1;
            dotind = dotind(end);
        else
            Filenames.DataFile = [file(1:dotind),'raw'];
            Filenames.ParameterFile = [file(1:dotind),'lab'];
        end
        
        try_path = which( 'MRecon.m' );
        try_path = try_path( 1:strfind( try_path, 'MRecon.m')-1);        
        if fopen( [try_path, 'recframe.exe'],'r') == -1
            cmd_str = which( 'recframe.exe' );
        else
            cmd_str = [try_path, 'recframe.exe'];
        end
       
        if isempty( cmd_str )
            error( 'recframe.exe not found. Please make sure that the location of recframe.exe is in the Matlab path');
        end        
        cmd_str = ['"',cmd_str, '" "', Filenames.DataFile, '" "', Filenames.ParameterFile, '" /D'];            
        system( cmd_str );               
        listfile = [Filenames.ParameterFile(1:dotind-1),'.list'];
        Parameter = listread(listfile);
        typ = unique(Parameter.Index.typ(:,2));
        for i = 1:length(typ)
            numtyp(i) = findstr(typ(i),t);
        end
        Parameter2read.typ = sort(numtyp);
        Parameter2read.mix = unique(Parameter.Index.mix);
        Parameter2read.dyn = unique(Parameter.Index.dyn);
        Parameter2read.card = unique(Parameter.Index.card);
        Parameter2read.echo = unique(Parameter.Index.echo);
        Parameter2read.loca = unique(Parameter.Index.loca);
        Parameter2read.chan = unique(Parameter.Index.chan);
        Parameter2read.extr1 = unique(Parameter.Index.extr1);
        Parameter2read.extr2 = unique(Parameter.Index.extr2);
        Parameter2read.ky = unique(Parameter.Index.ky);
        Parameter2read.kz = unique(Parameter.Index.kz);
        Parameter2read.aver = unique(Parameter.Index.aver);
        Parameter2read.rtop = unique(Parameter.Index.rtop);        
    otherwise
        Parameter2read = [];
        Parameter = [];
        DataFormat = [];
end
end

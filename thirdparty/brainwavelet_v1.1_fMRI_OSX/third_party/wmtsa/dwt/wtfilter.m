function [wtf] = wtfilter(wtfname, transform)
% wtfilter -- Define wavelet transform filter coefficients.
%
%****f* wmtsa.dwt/wtfilter
%
% NAME
%   wtfilter -- Define wavelet transform filter coefficients.
%
% SYNOPSIS
%  [wtf] = wtfilter(wtfname, transform)
%
% INPUTS
%   * wtfname    -- name of wavelet transform filter (string, case-insenstive).
%   * transform  -- name of wavelet transform  (string, case-insenstive).
%
% OUTPUTS
%   * wtf        -- wavelet tranform filter struct (wtf_s).
%
% SIDE EFFECTS
%
%
% DESCRIPTION
%   wtfilter returns a wtf_s struct with the wavelet (high-pass) and
%   scaling (low-pass) filter coefficients, and associated attributes. 
%
%   The wtf_s struct has fields:
%   * g         -- scaling (low-pass) filter coefficients (vector).
%   * h         -- wavelet (high-pass) filter coefficients (vector).
%   * L         -- filter length (= number of coefficients) (integer).
%   * Name      -- name of wavelet filter (character string).
%   * WTFclass  -- class of wavelet filters (character string).
%   * Transform -- name of transform (character string).
%
%   Typing wtfilter('list') displays a list of supported filters.
% 
%   Typing wtfilter('all') returns a struct array of wtf_s of all
%   supported filters.
%
%   The MODWT filter coefficients are calculated from the DWT filter 
%   coefficients:
%
%      ht = h / sqrt(2)
%      gt = g / sqrt(2)
%
%   The wavelet filter coefficients (h) are calculated from the scaling
%   filter coefficients via the QMF function (wmtsa_qmf).
% 
% USAGE
%  % Return a wtf struct.
%  wtf = wtfilter(wtfname, transform)   
%
%  % Display and return a list of available filter ames.
%  wtfnames = wtfilter('list')
%  
%  % Display and return a struct array of wtf_s for available filters.
%  wtfs = wtfilter('all', transform)
%
% ERRORS  
%   WMTSA:InvalidNumArguments
%   WMTSA:MissingRequiredArgument
%   WMTSA:InvalidArgumentValue
%
% EXAMPLE
%    wtf = wtfilter('LA8', 'modwt')
%    wtf = wtfilter('haar', 'dwt')
%
% ALGORITHM
%
%
% REFERENCES
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%     Time Series Analysis. Cambridge: Cambridge University Press.
%
% SEE ALSO
%  wtf_s, wtf_qmf
%
% TOOLBOX
%   wmtsa/dwt
%
% CATEGORY
%   Filters: Filters
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2005-02-28
%
% COPYRIGHT
%   (c) Charles R. Cornish 2005
%
% CREDITS
%   Based on the original function (myfilter.m) by Brandon Whitcher.
%
% REVISION
%   $Revision: 612 $
%
%***

% $Id: wtfilter.m 612 2005-10-28 21:42:24Z ccornish $

usage_str = ['Usage:  [wtf] = ', mfilename, ...
             ' (wtfname, transform)'];

ext_usage_str = [usage_str, '\n', ...
                 'Type ', mfilename, ...
                 '(''list'') for list of available filters.', '\n', ...
                 'Type ', mfilename, ...
                 '(''all'') for struct array of all available filters.', '\n' ...
                ];

%%   Check Input Arguments
error(nargerr(mfilename, nargin, [1:2], nargout, [0:1], ...
              1, sprintf(ext_usage_str), 'struct'));

%% Process 'list' option
if (strcmpi(wtfname, 'list'))
  wtfilters = available_wtfilters;
  wtfilter_names = {wtfilters.Name};
  disp('Wavelet transform filters are available for the following:');
  for ( i = 1:length(wtfilter_names))
    disp(wtfilter_names{i});
  end
  if (nargout > 0)
    wtf = wtfilter_names;
  end
  return
end

if (~exist('transform', 'var') || isempty(transform))
  error(['WMTSA:missingRequiredArgument'], ...
        encode_errmsg('WMTSA:missingRequiredArgument', wmtsa_err_table, 'transform'));
else
  if (isempty(strmatch(upper(transform), {'DWT', 'MODWT'}, 'exact')))
    error(['WMTSA:invalidArgumentValue'], ...
          encode_errmsg('WMTSA:invalidArgumentValue', wmtsa_err_table, 'transform', ...
                        transform));
  end
end

%% Process 'all' option
if (strcmpi(wtfname, 'all'))
  wtfilters = available_wtfilters;
  if (strcmpi(transform, 'MODWT'))
    for (i = 1:length(wtfilters))
      wtfilters(i).Transform = 'MODWT';
      wtfilters(i).g = wtfilters(i).g ./ sqrt(2);
      wtfilters(i).h = wtfilters(i).h ./ sqrt(2);
    end
  end
  wtf = wtfilters;
  return
end

hwtf = str2func(['wtf_', lower(wtfname)]);

try
  wtf = feval(hwtf);
catch
  error('WMTSA:invalidWaveletTransformFilter', ...
        ['Unrecognized or unsupported wavelet transform filter name ' , ...
         '(', wtfname, ').']);
end

if (strcmpi(transform, 'MODWT'))
  wtf.Transform = 'MODWT';
  wtf.g = wtf.g ./ sqrt(2);
  wtf.h = wtf.h ./ sqrt(2);
end

return

function wtfs = available_wtfilters
%   avaliable_filters -- Returns a struct array of the available filters.
  wtf_hlist = ...
  {
      @wtf_haar, ...
      @wtf_d4, ...
      @wtf_d6, ...
      @wtf_d8, ...
      @wtf_d10, ...
      @wtf_d12, ...
      @wtf_d14, ...
      @wtf_d16, ...
      @wtf_d18, ...
      @wtf_d20, ...
      @wtf_la8, ...
      @wtf_la10, ...
      @wtf_la12, ...
      @wtf_la14, ...
      @wtf_la16, ...
      @wtf_la18, ...
      @wtf_la20, ...
      @wtf_bl14, ...
      @wtf_bl18, ...
      @wtf_bl20, ...
      @wtf_c6, ...
      @wtf_c12, ...
      @wtf_c18, ...
      @wtf_c24
   };
  
  wtf = struct;
  for (i = 1:length(wtf_hlist))
    wtf = feval(wtf_hlist{i});
    wtfs(i) = wtf;
  end
  
return

% ******************************************************************************
% *** Extremal phase filters ***
% ******************************************************************************
function wtf = wtf_haar
  wtf.Name = 'Haar';
  wtf.g = [0.7071067811865475 0.7071067811865475];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 2;
  wtf.WTFClass = 'ExtremalPhase';
  wtf.Transform = 'DWT';
return
  
function wtf = wtf_d4
  wtf.Name = 'D4';
  wtf.g = [0.4829629131445341 0.8365163037378077 0.2241438680420134 ...
           -0.1294095225512603];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 4;
  wtf.WTFClass = 'ExtremalPhase';
  wtf.Transform = 'DWT';
  
function wtf = wtf_d6
  wtf.Name = 'D6';
  wtf.g = [0.3326705529500827 0.8068915093110928 0.4598775021184915 ...
           -0.1350110200102546 -0.0854412738820267 0.0352262918857096];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 6;
  wtf.WTFClass = 'ExtremalPhase';
  wtf.Transform = 'DWT';
return

function wtf = wtf_d8
  wtf.Name ='D8';
  wtf.g = [0.2303778133074431 0.7148465705484058 0.6308807679358788 ...
           -0.0279837694166834 -0.1870348117179132 0.0308413818353661 ...
           0.0328830116666778 -0.0105974017850021];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 8;
  wtf.WTFClass = 'ExtremalPhase';
  wtf.Transform = 'DWT';
return

function wtf = wtf_d10
  wtf.Name ='D10';
  wtf.g = [0.1601023979741930 0.6038292697971898 0.7243085284377729 ...
           0.1384281459013204 -0.2422948870663824 -0.0322448695846381 ...
           0.0775714938400459 -0.0062414902127983 -0.0125807519990820 ...
           0.0033357252854738];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 10;
  wtf.WTFClass = 'ExtremalPhase';
  wtf.Transform = 'DWT';
return

function wtf = wtf_d12
  wtf.Name ='D12';
  wtf.g = [0.1115407433501094 0.4946238903984530 0.7511339080210954 ...
           0.3152503517091980 -0.2262646939654399 -0.1297668675672624 ...
           0.0975016055873224 0.0275228655303053 -0.0315820393174862 ...
           0.0005538422011614 0.0047772575109455 -0.0010773010853085];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 12;
  wtf.WTFClass = 'ExtremalPhase';
  wtf.Transform = 'DWT';
return

function wtf = wtf_d14
  wtf.Name ='D14';
  wtf.g = [0.0778520540850081 0.3965393194819136 0.7291320908462368 ...
           0.4697822874052154 -0.1439060039285293 -0.2240361849938538 ...
           0.0713092192668312 0.0806126091510820 -0.0380299369350125 ...
           -0.0165745416306664 0.0125509985560993 0.0004295779729214 ...
           -0.0018016407040474 0.0003537137999745];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 14;
  wtf.WTFClass = 'ExtremalPhase';
  wtf.Transform = 'DWT';
return

function wtf = wtf_d16
  wtf.Name ='D16';
  wtf.g = [0.0544158422431049 0.3128715909143031 0.6756307362972904 ...
           0.5853546836541907 -0.0158291052563816 -0.2840155429615702 ...
           0.0004724845739124 0.1287474266204837 -0.0173693010018083 ...
           -0.0440882539307952 0.0139810279173995 0.0087460940474061 ...
           -0.0048703529934518 -0.0003917403733770 0.0006754494064506 ...
           -0.0001174767841248];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 16;
  wtf.WTFClass = 'ExtremalPhase';
  wtf.Transform = 'DWT';
return

function wtf = wtf_d18
  wtf.Name ='D18';
  wtf.g = [0.0380779473638791 0.2438346746125939 0.6048231236901156 ...
           0.6572880780512955 0.1331973858249927 -0.2932737832791761 ...
           -0.0968407832229524 0.1485407493381306 0.0307256814793395 ...
           -0.0676328290613302 0.0002509471148340 0.0223616621236805 ...
           -0.0047232047577520 -0.0042815036824636 0.0018476468830564 ...
           0.0002303857635232 -0.0002519631889427 0.0000393473203163];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 18;
  wtf.WTFClass = 'ExtremalPhase';
  wtf.Transform = 'DWT';
return

function wtf = wtf_d20
  wtf.Name ='D20';
  wtf.g = [0.0266700579005546 0.1881768000776863 0.5272011889317202 ...
           0.6884590394536250 0.2811723436606485 -0.2498464243272283 ...
           -0.1959462743773399 0.1273693403357890 0.0930573646035802 ...
           -0.0713941471663697 -0.0294575368218480 0.0332126740593703 ...
           0.0036065535669880 -0.0107331754833036 0.0013953517470692 ...
           0.0019924052951930 -0.0006858566949566 -0.0001164668551285 ...
           0.0000935886703202 -0.0000132642028945];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 20;
  wtf.WTFClass = 'ExtremalPhase';
  wtf.Transform = 'DWT';
return
  
% ******************************************************************************
% ***  Least asymmetric filters  ***
% ******************************************************************************
function wtf = wtf_la8
  wtf.Name ='LA8';
  wtf.g = [-0.0757657147893407 -0.0296355276459541 0.4976186676324578 ...
           0.8037387518052163 0.2978577956055422 -0.0992195435769354 ...
           -0.0126039672622612 0.0322231006040713];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 8;
  wtf.WTFClass = 'LeastAsymmetric';
  wtf.Transform = 'DWT';
return

function wtf = wtf_la10
  wtf.Name ='LA10';
  wtf.g = [0.0195388827353869 -0.0211018340249298 -0.1753280899081075 ...
           0.0166021057644243 0.6339789634569490 0.7234076904038076 ...
           0.1993975339769955 -0.0391342493025834 0.0295194909260734 ...
           0.0273330683451645];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 10;
  wtf.WTFClass = 'LeastAsymmetric';
  wtf.Transform = 'DWT';
return

function wtf = wtf_la12
  wtf.Name ='LA12';
  wtf.g = [0.0154041093273377 0.0034907120843304 -0.1179901111484105 ...
           -0.0483117425859981 0.4910559419276396 0.7876411410287941 ...
           0.3379294217282401 -0.0726375227866000 -0.0210602925126954 ...
           0.0447249017707482 0.0017677118643983 -0.0078007083247650];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 12;
  wtf.WTFClass = 'LeastAsymmetric';
  wtf.Transform = 'DWT';
return

function wtf = wtf_la14
  wtf.Name ='LA14';
  wtf.g = [0.0102681767084968 0.0040102448717033 -0.1078082377036168 ...
           -0.1400472404427030 0.2886296317509833 0.7677643170045710 ...
           0.5361019170907720 0.0174412550871099 -0.0495528349370410 ...
           0.0678926935015971 0.0305155131659062 -0.0126363034031526 ...
           -0.0010473848889657 0.0026818145681164];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 14;
  wtf.WTFClass = 'LeastAsymmetric';
  wtf.Transform = 'DWT';
return

function wtf = wtf_la16
  wtf.Name ='LA16';
  wtf.g = [-0.0033824159513594 -0.0005421323316355 0.0316950878103452 ...
           0.0076074873252848 -0.1432942383510542 -0.0612733590679088 ...
           0.4813596512592012 0.7771857516997478 0.3644418948359564 ...
           -0.0519458381078751 -0.0272190299168137 0.0491371796734768 ...
           0.0038087520140601 -0.0149522583367926 -0.0003029205145516 ...
           0.0018899503329007];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 16;
  wtf.WTFClass = 'LeastAsymmetric';
  wtf.Transform = 'DWT';
return

function wtf = wtf_la18
  wtf.Name ='LA18';
  wtf.g = [0.0010694900326538 -0.0004731544985879 -0.0102640640276849 ...
           0.0088592674935117 0.0620777893027638 -0.0182337707798257 ...
           -0.1915508312964873 0.0352724880359345 0.6173384491413523 ...
           0.7178970827642257 0.2387609146074182 -0.0545689584305765 ...
           0.0005834627463312 0.0302248788579895 -0.0115282102079848 ...
           -0.0132719677815332 0.0006197808890549 0.0014009155255716];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 18;
  wtf.WTFClass = 'LeastAsymmetric';
  wtf.Transform = 'DWT';
return

function wtf = wtf_la20
  wtf.Name ='LA20';
  wtf.g = [0.0007701598091030 0.0000956326707837 -0.0086412992759401 ...
           -0.0014653825833465 0.0459272392237649 0.0116098939129724 ...
           -0.1594942788575307 -0.0708805358108615 0.4716906668426588 ...
           0.7695100370143388 0.3838267612253823 -0.0355367403054689 ...
           -0.0319900568281631 0.0499949720791560 0.0057649120455518 ...
           -0.0203549398039460 -0.0008043589345370 0.0045931735836703 ...
           0.0000570360843390 -0.0004593294205481];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 20;
  wtf.WTFClass = 'LeastAsymmetric';
  wtf.Transform = 'DWT';
return

% ******************************************************************************
% ***  Best Localized filters  ***
% ******************************************************************************
function wtf = wtf_bl14
  wtf.Name ='BL14';
  wtf.g = [0.0120154192834842 0.0172133762994439 -0.0649080035533744 ...
           -0.0641312898189170 0.3602184608985549 0.7819215932965554 ...
           0.4836109156937821 -0.0568044768822707 -0.1010109208664125 ...
           0.0447423494687405 0.0204642075778225 -0.0181266051311065 ...
           -0.0032832978473081 0.0022918339541009];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 14;
  wtf.WTFClass = 'BestLocalized';
  wtf.Transform = 'DWT';
return 
  
function wtf = wtf_bl18
  wtf.Name ='BL18';
  wtf.g = [0.0002594576266544 -0.0006273974067728 -0.0019161070047557 ...
           0.0059845525181721 0.0040676562965785 -0.0295361433733604 ...
           -0.0002189514157348 0.0856124017265279 -0.0211480310688774 ...
           -0.1432929759396520 0.2337782900224977 0.7374707619933686 ...
           0.5926551374433956 0.0805670008868546 -0.1143343069619310 ...
           -0.0348460237698368 0.0139636362487191 0.0057746045512475];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 18;
  wtf.WTFClass = 'BestLocalized';
  wtf.Transform = 'DWT';
return 

function wtf = wtf_bl20
  wtf.Name ='BL20';
  wtf.g = [0.0008625782242896 0.0007154205305517 -0.0070567640909701 ...
           0.0005956827305406 0.0496861265075979 0.0262403647054251 ...
           -0.1215521061578162 -0.0150192395413644 0.5137098728334054 ...
           0.7669548365010849 0.3402160135110789 -0.0878787107378667 ...
           -0.0670899071680668 0.0338423550064691 -0.0008687519578684 ...
           -0.0230054612862905 -0.0011404297773324 0.0050716491945793 ...
           0.0003401492622332 -0.0004101159165852];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 20;
  wtf.WTFClass = 'BestLocalized';
  wtf.Transform = 'DWT';
return

% ******************************************************************************
% ***  Coiflet filters  ***
% ******************************************************************************
function wtf = wtf_c6
  wtf.Name ='C6';
  wtf.g = [-0.0156557285289848 -0.0727326213410511 0.3848648565381134 ...
           0.8525720416423900 0.3378976709511590 -0.0727322757411889];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 6;
  wtf.WTFClass = 'Coiflet';
  wtf.Transform = 'DWT';
return

function wtf = wtf_c12
  wtf.Name ='C12';
  wtf.g = [-0.0007205494453679 -0.0018232088707116 0.0056114348194211 ...
           0.0236801719464464 -0.0594344186467388 -0.0764885990786692 ...
           0.4170051844236707 0.8127236354493977 0.3861100668229939 ...
           -0.0673725547222826 -0.0414649367819558 0.0163873364635998];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 12;
  wtf.WTFClass = 'Coiflet';
  wtf.Transform = 'DWT';
return

function wtf = wtf_c18
  wtf.Name ='C18';
  wtf.g = [-0.0000345997728362 -0.0000709833031381 0.0004662169601129 ...
           0.0011175187708906 -0.0025745176887502 -0.0090079761366615 ...
           0.0158805448636158 0.0345550275730615 -0.0823019271068856 ...
           -0.0717998216193117 0.4284834763776168 0.7937772226256169 ...
           0.4051769024096150 -0.0611233900026726 -0.0657719112818552 ...
           0.0234526961418362 0.0077825964273254 -0.0037935128644910];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 18;
  wtf.WTFClass = 'Coiflet';
  wtf.Transform = 'DWT';
return

function wtf = wtf_c24
  wtf.Name ='C24';
  wtf.g = [-0.0000017849850031 -0.0000032596802369 0.0000312298758654 ...
           0.0000623390344610 -0.0002599745524878 -0.0005890207562444 ...
           0.0012665619292991 0.0037514361572790 -0.0056582866866115 ...
           -0.0152117315279485 0.0250822618448678 0.0393344271233433 ...
           -0.0962204420340021 -0.0666274742634348 0.4343860564915321 ...
           0.7822389309206135 0.4153084070304910 -0.0560773133167630 ...
           -0.0812666996808907 0.0266823001560570 0.0160689439647787 ...
           -0.0073461663276432 -0.0016294920126020 0.0008923136685824];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 24;
  wtf.WTFClass = 'Coiflet';
  wtf.Transform = 'DWT';
return

function wtf = wtf_c30
  wtf.Name ='C30';
  wtf.g = [-0.0000000951765727 -0.0000001674428858 0.0000020637618516 ...
           0.0000037346551755 -0.0000213150268122 -0.0000413404322768 ...
           0.0001405411497166 0.0003022595818445 -0.0006381313431115 ...
           -0.0016628637021860 0.0024333732129107 0.0067641854487565 ...
           -0.0091642311634348 -0.0197617789446276 0.0326835742705106 ...
           0.0412892087544753 -0.1055742087143175 -0.0620359639693546 ...
           0.4379916262173834 0.7742896037334738 0.4215662067346898 ...
           -0.0520431631816557 -0.0919200105692549 0.0281680289738655 ...
           0.0234081567882734 -0.0101311175209033 -0.0041593587818186 ...
           0.0021782363583355 0.0003585896879330 -0.0002120808398259];
  wtf.h = wmtsa_qmf(wtf.g);
  wtf.L = 30;
  wtf.WTFClass = 'Coiflet';
  wtf.Transform = 'DWT';
return


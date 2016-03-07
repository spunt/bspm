function lev = RestVolPath(Sub)
%%% Perform Standard Resting Pipline
%%% 1. Drop First Four Volumes
%%% 2. Perform Slice Time Correction
%%% 3. Perform Realignment
%%% 4. Perform Normalization to MNI Template
%%% 5. Smooth Data
%%% 6. Filter Data
%%% 7. Regress out motion and physiological variables
%%%
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%% Copyright (C) 2014,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.
% keyboard;
if exist(['res_ffBPS_motRes_ss_nn_st_dv_' Sub],'file')>0
    disp([Sub ' is all good']);
    lev = inf;
else
    if exist(['ffBPS_motRes_ss_nn_st_dv_' Sub],'file')>0
        lev = 7;
        RestingAna_ala_Buckner(['ffBPS_motRes_ss_nn_st_dv_' Sub],'fcParams',[8]);
        lev = inf;
    else
        if exist(['motRes_ff_ss_nn_st_dv_' Sub],'file')>0
            lev = 7;
            RestingAna_ala_Buckner(['motRes_ss_nn_st_dv_' Sub],'fcParams',[7 8]);
            lev = inf;
        else
            if exist(['ss_nn_st_dv_' Sub],'file')>0
                lev = 6;
                RestingAna_ala_Buckner(['ss_nn_st_dv_' Sub],'fcParams',[6 7 8]);
                lev = inf;
            else
                if exist(['nn_st_dv_' Sub],'file')>0
                    lev = 5;
                    RestingAna_ala_Buckner(['nn_st_dv_' Sub],'fcParams',[5 6 7 8]);
                    lev = inf;
                else
                    if exist(['realignment_params_st_dv_' Sub],'file')>0
                        lev = 4;
                        RestingAna_ala_Buckner(['st_dv_' Sub],'fcParams',[4 5 6 7 8]);
                        lev = inf;
                    else
                        if exist(['st_' Sub],'file')>0
                            lev = 3;
                            RestingAna_ala_Buckner(['st_dv_' Sub],'fcParams',[3 4 5 6 7 8]);
                            lev = inf;
                        else
                            if exist(['dv_' Sub],'file')>0
                                lev = 2;
                                RestingAna_ala_Buckner(['dv_' Sub],'fcParams',[2 3 4 5 6 7 8]);
                                lev = inf;
                            else
                                if exist(Sub,'file')>0
                                    lev = 1;
                                    RestingAna_ala_Buckner([Sub],'fcParams',[1 2 3 4 5 6 7 8]);
                                    lev = inf;
                                else
                                    lev = -1;
                                    disp(['No Nifti file in ' Sub '''s folder!'])
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

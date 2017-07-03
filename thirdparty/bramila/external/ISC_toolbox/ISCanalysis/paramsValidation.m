function [Params,flag] = paramsValidation(Params)

h.Pub = Params.PublicParams;
h.checkboxFreq = true;
h = validateDataAndParams(h);
disp(' ')
if h.validFlag
    disp('Type "runAnalysis(Params)" to start analysis.')
    Pr.PublicParams = h.Pub;
    Pr.PrivateParams = h.Priv;
    Params = Pr;
else
    disp('You must set parameters correctly before running the analysis.')
end


if nargout == 2
    flag = h.validFlag;
end
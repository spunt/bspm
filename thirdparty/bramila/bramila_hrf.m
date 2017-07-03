function hrfOut = bramila_hrf(par1,par2,par3,par4)
%Produces a canonical two Gamma hemodynamic response function (HRF)
%
%Takes four parameters:
%par1 time step length in seconds
%par2 weight for the second Gamma (e.g. 0 for just one gamma funtion)
%par3 shape parameter for first Gamma probability density function (pdf)
%par4 shape parameter for second Gamma pdf
%
%defaults:
%par1 = 2
%par2 = 1
%par3 = 6
%par4 = 16
%
%The scale parameter for both Gamma pdfs is 1/par3
%length of the output HRF is 30 seconds

if nargin<1 || isempty(par1)
    par1=2;
end;
if nargin<2 || isempty(par2)
    par2=1;
end;
if nargin<3 || isempty(par3)
    par3=6;
end;
if nargin<4 || isempty(par4)
    par4=16;
end;


timepoints=0:par1:30-par1;

hrfOut=gampdf(timepoints,par3,1)-par2*gampdf(timepoints,par4,1)/par3;

end

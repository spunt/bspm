function[z]=fisherz(r)
%FISHERZ Fisher's Z-transform.
%   Z = FISHERZ(R) returns the Fisher's Z-transform of the correlation
%   coefficient R.

%20080103, Thomas Zoeller (tzo@gmx.de)

r=r(:);
z=.5.*log((1+r)./(1-r));
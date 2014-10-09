function cofs=d2boxcof
% cofs=D2BOXCOF
%
% Puts all the required coefficients for the 2-tap Daubechies wavelet
% transform on the interval in a structure array.
%
% SEE ALSO: D2BOXSTEP, D2BOXSTEPI
%
% Last modified by fjsimons-at-alum.mit.edu, 08/24/2010

% Get the product filter coefficients
[h0,f0]=wc('Daubechies',1);
% Get the lowpass (H0) and highpass (H1) coefficients of the ANALYSIS
% Get the lowpass (F0) and highpass (F1) coefficients of the SYNTHESIS
[H0,H1,F0,F1]=prodco(f0,h0);

% Collect all of these filter coefficients into a structure
cofs.H0=H0; cofs.H1=H1; cofs.F0=F0; cofs.F1=F1; 

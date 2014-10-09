function [rad,rho,PV,PH,SV,SH,QK,QMu,eta]=getmodel(mod1,dir1,catflag)
% [rad,rho,PV,PH,SV,SH,QK,QMu,eta]=GETMODEL(mod1,dir1,catflag)
%
% Loads an Earth model, suitable for input to normal-mode codes.
%
% INPUT:
%
% mod1        The model file name [default: 'prem800']
% dir1        The directory name [default: $IFILES/EARTHMODELS/]
% catflag     0 YANNOS format by Yann Capdeville
%             1 MINEOS format by Masters, Woodhouse & Gilbert
%
% OUTPUT:
%
% rad         Radius [m]
% rho         Density [kg/m3]
% PV          Vertically-polarized P-wave speed [m/s]
% PH          Horizontall- polarized P-wave speed [m/s]
% SV          Vertically-polarized P-wave speed [m/s]
% SH          Horizontally-polarized P-wave speed [m/s]
% QK          Compressional quality factor
% QMU         Shear quality factor
% eta         Anisotropic parameter
%
% EXAMPLES: See the defaults inside this function
%
% SEE ALSO:
%
% GETSPHEROIDAL, GETTOROIDAL, EQPOTENTIAL
%
% Last modified by jchawtho-at-princeton.edu, 07/01/2007
% Last modified by efwelch-at-princeton.edu, 07/21/2010
% Last modified by fjsimons-at-alum.mit.edu, 07/10/2012

defval('catflag',1);

% Reads the header
switch catflag
 case 0
  % YANNOS model
  % An isotropic prem
  defval('mod1','prem800')
  defval('dir1',fullfile(getenv('IFILES'),'EARTHMODELS','MODES','YANNOS'))
  
  % Model name
  mname=fullfile(dir1,sprintf('%s',mod1));
  
  mfid=fopen(mname);
  
  C=textscan(mfid,'%s',1);
  modname=C{1};
  
  C=textscan(mfid,'%d%f%d',1);
  layered=C{1}; reffreq=C{2}; isanisotropic=C{3};
  
  C=textscan(mfid,'%d%d%d%d',1);
  nlay=C{1}; inncorelay=C{2}; cmblay=C{3}; numoceanlay=C{4};
 case 1 
  % MINEOS model (the input file for ./minos_bran)
  % An anisotropic prem
  defval('mod1','aniprem808')
  defval('dir1',fullfile(getenv('IFILES'),'EARTHMODELS','MODES','MINEOS-1.0.0'))
  % An anisotropic prem with the ocean removed
  defval('mod1','noocean_aniprem808')
  defval('dir1',fullfile(getenv('IFILES'),'EARTHMODELS','MODES','MINEOS-1.0.2'))
  
  % Model name
  mname=fullfile(dir1,sprintf('%s.dk',mod1));
  
  mfid=fopen(mname);
  
  modname=fgetl(mfid);
  
  C=textscan(mfid,'%d%f%d',1);
  format=C{1}; reffreq=C{2}; isanisotropic=C{3};
  
  C=textscan(mfid,'%d%d%d',1);
  nlay=C{1}; icblay=C{2}; cmblay=C{3};
end

% Read the actual model
C=textscan(mfid,'%f%f%f%f%f%f%f%f%f');

% Radius, m
rad=C{1};

% Density, kg m^-3
rho=C{2};
 
if nargout>2
  % Velocities, m/s
  PV=C{3};
  SV=C{4};
  QK=C{5};
  QMu=C{6};
  PH=C{7};
  SH=C{8};
  eta=C{9};
end

fclose(mfid);

disp(sprintf('Using model %s',mname))


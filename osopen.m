function [fid0,fid1,fid2,fid3,fmt1,fmt2,fmt3,fmtf,fmte,fmtd,fmtc,fmtb,fmta]=...
    osopen(np,npp)
% [fid0,fid1,fid2,fid3,fmt1,fmt2,fmt3,fmtf,fmte,fmtd,fmtc,fmtb,fmta]=...
%     OSOPEN(np,npp)
%
% Opens a bunch of diagnostic files and returns a ton of format strings
%
% INPUT:
%
% np     The number of parameters to solve for (e.g. 3, 5 or 6)
% npp    The number of unique entries in an np*np symmetric matrix
%
% SEE ALSO:
%
% OSLOAD, DIAGNOS (with which it needs to match!)
%
% Last modified by fjsimons-at-alum.mit.edu, 10/02/2014

% Who called? Work this into the filenames
[~,n]=star69;

% Ouput files, in parallel might be a jumble
% The thruth and the theoretical covariances
fid0=fopen(sprintf('%s_thzro_%s',n,date),'w');
% The estimates
fid1=fopen(sprintf('%s_thhat_%s',n,date),'a+');
% The initial guesses
fid2=fopen(sprintf('%s_thini_%s',n,date),'a+');
% The collected optimization diagnostics, replicated some above
fid3=fopen(sprintf('%s_diagn_%s',n,date),'a+');

% Output formatting for the estimation parameters
if np==3
  %                           s2    nu    rho
  fmt1=[                   '%9.3e %6.3f %6.0f\n'];
elseif np==5
  %                D    f2    s2    nu    rho
  fmt1=[      '%12.6e %6.3f %9.3e %6.3f %6.0f\n'];
elseif np==6
  %         D    f2      r    s2    nu    rho
  fmt1=['%12.6e %6.3f %6.3f %9.3e %6.3f %6.0f\n'];
end

% Output formatting for the simulation parameters
if np>=5
  %     DEL     g  z2  dydx  NyNx  blurs
  fmt2='%i %i %5.2f %i %i %i %i %i %i\n';
else
  %     dydx  NyNx  blurs
  fmt2='%i %i %i %i %i\n';
end

% For the time, exit flag, iterations 
fmta='%3i %3i %3i\n';
% For the likelihood and first-order optimality
fmtb='%15.8e %15.8e\n';
% For the scale
fmtc=[repmat('%15.0e ',1,np) '\n'];
% For the score
fmtd=[repmat('%15.8e ',1,np) '\n']; 
% For the Hessian, with npp unique elements:
fmte=repmat([repmat('%15.12f ',1,npp/3) '\n'],1,3); 
% For the unscaled covariance matrix
fmtf=repmat([repmat('%19.12e ',1,npp/3) '\n'],1,3);

% Lumps some of the formats together
fmt3=[fmta fmtb fmtc fmtd fmte];

function [fid0,fid1,fid2,fid3,fmt1,fmt2,fmt3,fmtf,fmte,fmtd,fmtc,fmtb,fmta]=...
    osopen(np)
% [fid0,fid1,fid2,fid3,fmt1,fmt2,fmt3,fmtf,fmte,fmtd,fmtc,fmtb,fmta]=...
%     OSOPEN(np,npp)
%
% Opens a bunch of diagnostic files and returns a ton of format strings
%
% INPUT:
%
% np     The number of parameters to solve for (e.g. 3, 5 or 6)
%
% SEE ALSO:
%
% OSLOAD, DIAGNOS (with which it needs to match!)
%
% Last modified by fjsimons-at-alum.mit.edu, 10/22/2014

% Who called? Work this into the filenames
[~,n]=star69;

% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;

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
  %     DEL     g  z2  dydx  NyNx  blurs kiso
  fmt2='%i %i %5.2f %i %i %i %i %i %i %f\n';
else
  %     dydx  NyNx  blurs kiso
  fmt2='%i %i %i %i %i %f\n';
end

% For the time, exit flag, iterations 
fmta='%3i %3i %3i\n';
% For the likelihood, first-order optimality, and moments
fmtb='%15.8e %15.8e %15.8e %15.8e %15.8e\n';
% For the scale
fmtc=[repmat('%15.0e ',1,np) '\n'];
% For the score, the gradient of the misfit function
fmtd=[repmat('%15.8e ',1,np) '\n']; 
% For the Hessian of the misfit function, with npp unique elements 
fmte=repmat([repmat('%15.12f ',1,npp/3) '\n'],1,3);
% For the unscaled Hessian-derived covariance matrix, or 
% for the unscaled theoretical covariance matrix
fmtf=repmat([repmat('%19.12e ',1,npp/3) '\n'],1,3);

% Lumps some of the formats together
fmt3=[fmta fmtb fmtc fmtd fmte fmtf];

function varargout=jxm2plm(fname)
% lmcosi=jxm2plm(fname)
% [hlm,dels,dems]=jxm2plm(fname)
%
% Converts a spherical harmonics file from Jerry X. Mitrovica 
%
% INPUT:
%
% fname     A complete file name string
%
% OUTPUT:
%
% lmcosi    A matrix of expansion coefficients in a format suitable for
%           plotting and analysis by PLOTPLM, PLM2XYZ, etc., OR
%
% hlm       The coefficients in a linear vector [indeks(lmcosi(:,3:4)',':')]
% dels      The spherical harmonic degrees [this would be lmcosi(:,1)]
% dems      The spherical harmonic orders [this would be lmcosi(:,2)]
%
% Last modified by fjsimons-at-alum.mit.edu, 07/23/2010

% Load the file
fid=fopen(fname,'r');
% Read the header
head=fgetl(fid);
LT=sscanf(head,'%f',2);
% Spherical harmonic bandwidth
L=LT(1);
% Time
T=LT(2);
% Read the coefficients
COF=fscanf(fid,'%f');
fclose(fid);
% Check that the length is as I expect it
difer(length(COF)-2*addmup(L)-2*(L/2+1),[],[],NaN)

% Now rearrange
mods=(L-[0:L]+1);
% The indices of the blank pairs
modu=cumsum(mods+mod(mods,2));
modb=modu(~~mod(mods,2));
blanx=sort([2*modb-1 2*modb]);
% Check these are indeed blanks
difer(COF(blanx),[],[],NaN)
% Remaining are the non-blanks
COF=skip(COF,blanx);
% Which should be the same number as expected
difer(length(COF)-2*addmup(L),[],[],NaN)

% Now split into the real and imaginary parts
C=COF(1:2:end);
S=COF(2:2:end);
% The first L+1 should be zero
difer(S(1:L+1),[],[],NaN)

% Create a blank array with FJS format
[dems,dels,mz,lmcosi]=addmon(L);
% Fake arrival to align with readshds
hlm=lmcosi(:,3:4);
modm=[cumsum([0 mods(1:end-1)])+1 addmup(L)+1];
% Populate this with the JXM coefficients
for m=0:L
  % Stick in the "cosine" coefficients
  hlm(addmup(m-1:L-1)+m+1,1)=C(modm(m+1):modm(m+2)-1);
  % Stick in the "sine" coefficients
  hlm(addmup(m-1:L-1)+m+1,2)=S(modm(m+1):modm(m+2)-1);
end

% Some conventions need to be realigned
CSC=(-1).^dems;
dom=sqrt(2-(dems==0));
hlm(:,1)=hlm(:,1).*CSC.*dom;
hlm(:,2)=hlm(:,2).*CSC.*dom;

% Preserve the natural order of things
lmcosi(:,3:4)=hlm;

% Fake arrival to align with readshds
hlm=hlm';
hlm=hlm(:);

if nargout==1
  varargout={lmcosi};
end
if nargout==3
  varargout{1}=hlm;
  varargout{2}=dels;
  varargout{3}=dems;
end

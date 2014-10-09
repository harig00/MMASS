function [cumnum,gcd,azm,celnr,c11,cmn,lat,...
	  dlon,nmr,lonlat,indi,refarea,basdep,basfun]=...
    lcsg(diro,genx)
% [cumnum,gcd,azm,celnr,c11,cmn,lat,...
%   dlon,nmr,lonlat,indi,refarea,basdep,basfun]=...
%      LCSG(diro,genx)
% 
% In a raytracing directory 'diro', opens up the five files
% created by 'surftrace' and returns their content.
% 'L.bin' Cumulative nr of entries per row (event)    => cumnum
% 'C.bin' Cell (column) number of each entry          => celnr
% 'S1.bin' Values of all entries (path lengths in km) => gcd
% 'S2.bin' Values of all entries (azimuths)           => azm
% 'lonlat'                                            => lonlat
%
% Also returns a matrix 'indi' with, for every path, the
% begin and end indices of the long vectors to be considered.
% Also returns the basis functions from the 'inamak3d' file.
%
% Gets gridinfo (does gridlod2) from $GRIDS
% 'gridinfo'                                          => c11,cmn,lat,dlon,nmr
%
% SEE ALSO: LCS2SPR, LCSD
%
% Last modified by fjsimons-at-alum.mit.edu, 08/23/2007

fid1=fopen(fullfile(diro,'L.bin'),'r'); 
fid2=fopen(fullfile(diro,'S1.bin'),'r');
fid3=fopen(fullfile(diro,'S2.bin'),'r');
fid4=fopen(fullfile(diro,'C.bin'),'r'); 

if any([fid1 fid2 fid3 fid4]==[-1 -1 -1 -1])
  error('Files not found')
end

[c11,cmn,lat,dlon,nmr,lej,dieptes,depth,...
 disc,basdep,basfun,X,Y,XI,YI,refarea]=...
    gridlod2(genx);

load(fullfile(diro,'lonlat'))

cumnum=fread(fid1,'int32','b');
gcd=fread(fid2,'float32','b');
azm=fread(fid3,'float32','b');
celnr=fread(fid4,'int32','b');

indi=[0 ; cumnum];

indi=[indi(1:end-1)+1 indi(2:end)];

fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);

[m1,n1]=size(cumnum);
[m2,n2]=size(indi);
[m3,n3]=size(lonlat);
if any(diff([m1 m2 m3])) ; error ; end
[m1,n1]=size(azm);
[m2,n2]=size(celnr);
[m3,n3]=size(gcd);
if any(diff([m1 m2 m3])) ; error ; end



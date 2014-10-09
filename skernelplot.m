function skernelplot(kdiro,dshell,ftype)
% SKERNELPLOT(kdiro,dshell,ftype)
%
% Plots "banana-dohnut" sensitivity kernels on the cubed sphere.
%
% INPUT:
%
% kdiro      A directory name with a kernel in it
% dshell     An integer identifying the layer number
% ftype      1 Identifies the first ever metafile used
%            2 Identifies the seconnd metafile type used
%
% EXAMPLE:
%
% Some tests of the directories with Stephen's kernels at three depths
%
% diro='/u/fjsimons/MyPapers/WAVELETS/StephenJudd';
% kdiros=ls2cell(fullfile(diro,'kernelsEtc'));
% for index=1:length(kdiros)
%  skernelplot(fullfile(diro,'kernelsEtc',kdiros{index}),100)
%  skernelplot(fullfile(diro,'kernelsEtc',kdiros{index}),110)
%  skernelplot(fullfile(diro,'kernelsEtc',kdiros{index}),128)
% end
%
% Last modified by fjsimons-at-alum.mit.edu, 11/20/2009

% Supply defaults
defval('kdiro',...
       '/u/fjsimons/MyPapers/WAVELETS/StephenJudd/kernelsEtc/000220Cuba.krnl')
defval('dshell',100)
defval('filetype',2);

% Read the kernel into memory
m=readKernelMex(kdiro);

% Read the depth divisions of the cubed sphere
load(fullfile(getenv('IFILES'),'EARTHMODELS','RITSEMA','xpcube.mat'))
% The result is a variable 'xp' with 129 depths in the first column
[d,i]=sort(xp(:,1),1,'descend'); 
% Now you get the depths, i.e. 2889.6 2867.0 2844.4 etc. up to 0 
depths=xp(i,1);

% Load the metadata also
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch filetype
  case 1
   fid=fopen(fullfile(kdiro,'meta.data'));
   fgetl(fid); fgetl(fid);
   [ss,slat,slon,sdep]=strread(fgetl(fid),'%s %f %f %f');
   [rs,rlat,rlon,rs,rs]=strread(fgetl(fid),'%s %f %f %s %f');
   rdep=0; fgetl(fid); fgetl(fid); 
   [cg,cg,nchunk,nxi,neta,ndep]=strread(fgetl(fid),'%s %s %f %f %f %f');
   fclose(fid);
 case 2
   fid=fopen(fullfile(kdiro,'meta.data'));
   fgetl(fid); fgetl(fid);
   [ss,slat,slon,sdep]=strread(fgetl(fid),'%s %f %f %f');
   [rs,rlat,rlon,rdep]=strread(fgetl(fid),'%s %f %f %f');
   fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); 
   fgetl(fid); fgetl(fid); 
   [cg,cg,nchunk,nxi,neta,ndep]=strread(fgetl(fid),'%s %s %f %f %f %f');
   fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that the metdata corresponds to the model read
difer([ndep nxi neta nchunk]-size(m),[],[],NaN)
difer(ndep-length(depths),[],[],NaN)

% Make the actual figure
clf
plotoncube(squeeze(double(m(dshell,:,:,:))),'2D'); xy=axis;
plotcont([],[],9); xlim(xy(1:2)+[-1 1]*range(xy(1:2)/100))
hold on
[xic,etac]=sphere2cube(slon,slat);
[spc,pgc]=plotonchunk(xic,etac);  delete(cat(1,pgc{:})) 
hold on
[xic,etac]=sphere2cube(rlon,rlat);
[rpc,pgc]=plotonchunk(xic,etac);  delete(cat(1,pgc{:})) 
hold off
set([spc(~isnan(spc)) rpc(~isnan(rpc))],'MarkerS',10)

% Cosmetics and diagnostics
kname=pref(suf(kdiro,'/'));
title(kname)
xloc=-0.75;
t(1)=text(xloc,3.5,...
    sprintf('lon = %4.1f%s ; lat = %4.1f%s ; dep = %4.1f km',...
	    slon,str2mat(176),slat,str2mat(176),sdep));
t(2)=text(xloc,3.25,...
    sprintf('lon = %4.1f%s ; lat = %4.1f%s ; dep = %4.1f km',...
	    rlon,str2mat(176),rlat,str2mat(176),rdep));

t(3)=text(xloc,3.0,...
    sprintf('rendered depth shell = %5.1f km',depths(dshell)));

% Cosmetics and print suggestion
fig2print(gcf,'landscape')
figdisp([],sprintf('%s_%i',kname,dshell),[],1)


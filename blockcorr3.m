function varargout=blockcorr3(diro,outp,strt)
% [sr,sc,corm,m,n,a]=blockcorr3(diro,outp,strt)
%
% Figures out the shifts necessary to merge a series of images.
%
% INPUT:
%
% diro       Directory in which the [*bmp] files reside
% outp       Output (create or append) filename string [default: blox]
% strt       The number of the first image to be considered [default: 1]
% 
% OUTPUT:
%
% sr,sc      Row/col shifts (plus is down and right) of im2's onto im1's
% corm       Correlation coefficients at the maximized location
% m,n        Row/col dimensions of the series of images
% a          A cell array with the file names used in this sequence
%
% EXAMPLE:
%
% blockcorr3('/home/fjsimons/STUDENTS/CatherineRose',[],153)
% blockcorr3('/home/fjsimons/STUDENTS/CatherineRose',[],248)
%
% SEE ALSO:
%
% This is the workhorse, and uses or leads to the following:
% BLOCKCORR, BLOCKCORR2, INNERMERGE
%
% Last modified by fjsimons-at-alum.mit.edu, 06/10/2008

% Default directory
defval('diro','/home/fjsimons/STUDENTS/CatherineRose');
defval('outp','/home/fjsimons/STUDENTS/CatherineRose/blox');
defval('strt',1);

% First, on the regular images - after you're done, could check that the
% remaining correlations are all 1 and point to zero shifts, or
% thereabouts. This verification can be done in a very small extent 

% Get all the file names - in bmp [and/or jpg]
% a=ls2cell(fullfile(diro,'s*bmp'));
a=ls2cell(fullfile(diro,'scn*'));

% In this code we're only loading one image at the time, to allow for
% unequally sized but same-resolution ones.

% Get a first image, convert to grey, initialize variables
im2=greyimg(getimg(diro,a,strt));
[m(1),n(1)]=size(im2);
sr(1)=0; sc(1)=0; corm(1)=0;

% Open a file and start doing the analysis
fid=fopen(outp,'a'); % When append, watch the first file
if strt==1
  fwrite(fid,sprintf('%s %5i %5i %3i %3i %8.4f\n',a{1},m(1),n(1),0,0,1));
end
for in=strt:length(a)-1
  % Get first image from stack
  im1=im2;
  % Get second image, convert to grey
  im2=greyimg(getimg(diro,a,in+1));
  [m(in+1),n(in+1)]=size(im2);

  % Now do the block cross-correlation - 0 at the end means no plots
  [sr(in+1),sc(in+1),corm(in+1)]=blockcorr2(im1,im2,1);

  colormap(gray(256))
  drawnow
  
  % Write output on the fly since this is a long process
  % img=str2num(strtok(pref(a{in+1}),'scn'));
  fwrite(fid,sprintf('%s %5i %5i %3i %3i %8.4f\n',a{in+1},m(in+1),n(in+1),...
		     sr(in+1),sc(in+1),corm(in+1)));
end
fclose(fid);

% Ouput
vars={sr,sc,corm,m,n,a};
varargout=vars(1:nargout);

% Later reinvent these images by shifting them! This happens in INNERMERGE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts an image to grey scale
function img=greyimg(im)
[m,n,o]=size(im);

% If not initializing, must convert to double
reds=nan(m,n);
greens=nan(m,n);
blues=nan(m,n);
greys=nan(m,n);

% Isolate channels - explicit indexing if not -> double
reds(:,:)=reshape(indeks(im,1:m*n),m,n);
greens(:,:)=reshape(indeks(im,[1:m*n]+m*n),m,n);
blues(:,:)=reshape(indeks(im,[1:m*n]+2*m*n),m,n);

% Convert to grey scale
img=(reds+blues+greens)/3;

% Demo function - get a real image
function im=getimg(diro,aa,ins)
% Get an image
getthis=fullfile(diro,aa{ins});
im=imread(getthis);
display(sprintf('Loading %s',getthis))

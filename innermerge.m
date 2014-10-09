function innermerge(sr,sc,corm,m,n,a)
% INNERMERGE(sr,sc,corm,m,n,a)
%
% Calculates the largest common overlap region for a series of images
% block-correlated using BLOCKCORR3 and writes new images
%
% INPUT:
%
% sr         A string with a filename containing the following, OR 
% sr,sc      Row/col shifts (plus is down and right) of im2's onto im1's
% corm       Correlation coefficients at the maximized location
% m,n        Row/col dimensions of the series of images
% a          A cell array with the file names used in this sequence
%
% EXAMPLE:
%
% innermerge
%
%
% Last modified by fjsimons-at-alum.mit.edu, 05/08/2008

defval('sr','blox');
defval('diro','/home/fjsimons/STUDENTS/CatherineRose');

if isstr(sr)
  % Anything else (fread, fscanf, etc) won't work - I tried
  [a,m,n,sr,sc,corm]=textread(fullfile(diro,sr),...
			      '%s %n %n %n %n %n\n');
  % Make consistent with what's in the directory
end

% Check you've got all of them
difer(2*length(a)-length(sr)-length(sc))

% First of all, reference all the offsets to the first image
% These are the starting points of all the images in the coordinate
% system of the first
sr=cumsum(sr(:));
sc=cumsum(sc(:));

% Figure out the maximum inner dimensions with respect to the first
% These are the end points of all the images in the coordinate system of
% the first
er=sr(:)+m(:);
ec=sc(:)+n(:);

% This is the range of rows and columns we want of the first image, in
% the coordinates of the first image
rowr=max(1+sr):min(er);
colr=max(1+sc):min(ec);

% For the next few images
for in=1:length(a)
  % Read the images, one by one
  im=imread(fullfile(diro,a{in}));
  % And now the range that we want is the same as the first image, except
  % in the new coordinates, and making sure we don't exceed the allowable
  % values 
  imp=im(intersect(rowr-sr(in),1:m(in)),...
	 intersect(colr-sc(in),1:n(in)),:);
  % Check the consistent size
  % disp(sprintf('Making a new image at %i%s%i pixels',...
  % size(imp,1),'x',size(imp,2)))
  % Resave the images
  imwrite(imp,fullfile(diro,sprintf('n%s',a{in})))
end



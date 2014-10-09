function varargout=blockcorr2(im1,im2,plots)
% [sr,sc,coma]=blockcorr2(im1,im2,plots)
%
% Finds the shifts required of im2 to be aligned with im1
% Both im1 and im2 are very large; the algorithm works with a portion of
% the first image and shifts it about the second image. Work on
% single-channel colors only, e.g. the grey values.
%
% INPUT:
%
% im1        The first single-channel image [mxnx1]
% im2        The second image, of exactly the same size
% plots      0 Only produce numerical output [default]
%            1 Make diagnostic plots on the way
%
% OUTPUT:
%
% sr,sc      Row/col shift (plus is down and right) to im2 to align with im1
% coma       Maximum of the cross-correlation function
%
% EXAMPLE:
%
% blockcorr2('demo1') % For two different images
% blockcorr2('demo2') % To make zero shifts are also detected
%
% Last modified by fjsimons-at-alum.mit.edu, May 1st, 2008

if ~isstr(im1) % Not the demo
  defval('plots',0)
  % Check compatibility of sizes of both single-channel images
  % difer(size(im1)-size(im2))
  [m,n]=size(im1);

  % The row/col index of the center(s) of the tile(s) in im1
  cr=round(m/2);
  cc=round(n/2);
  % The row and column extent(s) (size[s]) of the tile(s) in im1 
  er=round(m/10);
  ec=round(n/10);
  % The row/col index of the centers of the correlation region in im2
  % Best to start from the same locations
  cor=cr;
  coc=cc;
  % The row/col extent of the correlation regions in im2
  % Best to make this a small value based on the subtile
  % Experiment with this! This is going to need to vary per experiment.
  cer=round(er/1.5);
  cec=round(ec/1.5);

  % Perform the measurement as many times as you specified
  for in=1:length(cr)
    % Now perform the actual correlation with these values
    [sr(in),sc(in),corm(in),corv,maxi,maxr,maxc,sub1,sub2]=...
	blockcorr(im1,im2,[cr(in) cc(in)],[er(in) ec(in)],...
		  [cor(in) coc(in)],[cer(in) cec(in)]);
    
    [ri,ci]=ind2sub(fliplr(size(corv)),maxi);
    % The below is correct despite the appearances        
    if ri==1 || ci==1 || ri==size(corv,2) || ci==size(corv,1) || ...
	  corm(in)<0.5 % What if it's just bad?
      error('Widen the cross-correlation search')
    end
    
    if plots==1
      % As part of this loop maybe you want to plot the subimages which
      % should look perfectly aligned (and thus the big image will)
      clf
      [ah,ha]=krijetem(subnum(1,3));
      % Subimage in im1
      axes(ah(1))
      imagesc(sub1); axis image; longticks
      set(ah(1),'xtick',[1 ec(in)])
      set(ah(1),'ytick',[1 er(in)])
      title(sprintf('%i%s%i tile of im1 on (%i,%i)',...
	    er(in),'\times',ec(in),cr(in),cc(in)))
      % Subimage in im2
      axes(ah(2))
      imagesc(sub2); axis image; longticks
      set(ah(2),'xtick',[1 ec(in)])
      set(ah(2),'ytick',[1 er(in)])
      title(sprintf('%i%s%i tile of im2 on (%i,%i)',...
	    er(in),'\times',ec(in),maxr,maxc))
      
      % Correlation values - the maximum should be well defined            
      axes(ah(3))
      % Plot the negative to get black where the maximum is
      imagesc(-corv); axis image; longticks
      set(ah(3),'xtick',[1 cec(in)+~mod(cec(in),2)])
      set(ah(3),'ytick',[1 cer(in)+~mod(cer(in),2)])
      title(sprintf('%i%s%i cross-correlation',...
		    cer(in)+~mod(cer(in),2),'\times',...
		    cec(in)+~mod(cec(in),2)))
      xlabel(sprintf('Shift im2 by %+i rows and %+i cols',sr,sc))
      hold on

      % The notation isn't optimal but the point is it works
      s=plot(ri,ci,'o'); set(s,'MarkerFaceC','y','MarkerEdgeC','r')
      if length(cr)>1
	pause
      end
      display(sprintf('Shift im2 by %+i rows and %+i cols',sr,sc))
    end
  end
  % Now look at the results - these should be pretty consistent
  % and make an aggregate... or maybe just stick to one thing

  % Output
  vars={sr,sc,corm};
  varargout=vars(1:nargout);
  
elseif strmatch(im1,'demo1')
  % Get two test images in full color
  [im1,im2]=getimg;
  % Convert to a single channel, here the grey values
  [im1g,im2g]=greyimg(im1,im2);
  % Perform the block correlation on some portion(s)
  [mf,nf,o]=size(im1);
  [m,n]=size(im1g);
  % Check it's all well
  difer(m-mf)
  difer(n-nf)
  % Perform the operation which itself calls BLOCKCORR - with plots
  [sr,sc]=blockcorr2(im1g,im2g,1);
elseif strmatch(im1,'demo2')
  % Get two test images in full color
  im1=getimg; im2=im1;
  % Convert to a single channel, here the grey values
  [im1g,im2g]=greyimg(im1,im2);
  % Perform the block correlation on some portion(s)
  [mf,nf,o]=size(im1);
  [m,n]=size(im1g);
  % Check it's all well
  difer(m-mf)
  difer(n-nf)
  % Perform the operation which itself calls blockcorr - with plots
  [sr,sc]=blockcorr2(im1g,im2g,1);
  % Should be both zero; the images were identical!
  difer(sr+sc)
end

% Converts the given images to grey scales and base algorithms on it
function [im1g,im2g]=greyimg(im1,im2)
[m,n,o]=size(im1);
reds=nan(m,n,2);
greens=nan(m,n,2);
blues=nan(m,n,2);
greys=nan(m,n,2);
redit={im1,im2};
for in=1:2
  reds(:,:,in)=reshape(indeks(redit{in},1:m*n),m,n);
  greens(:,:,in)=reshape(indeks(redit{in},[1:m*n]+m*n),m,n);
  blues(:,:,in)=reshape(indeks(redit{in},[1:m*n]+2*m*n),m,n);
  greys(:,:,in)=(reds(:,:,in)+blues(:,:,in)+greens(:,:,in))/3;
end
% Return the two gray scale images
im1g=greys(:,:,1);
im2g=greys(:,:,2);

% Demo function - get two real images
function [im1,im2]=getimg
% Get two images (or one)
diro='/home/fjsimons/STUDENTS/CatherineRose';
a=ls2cell(fullfile(diro,'*bmp'));
im1=imread(fullfile(diro,a{1}));
if nargout>1
  im2=imread(fullfile(diro,a{2}));
end

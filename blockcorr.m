function varargout=blockcorr(im1,im2,crcc,erec,corcoc,cercec)
% [sr,sc,corm,corv,maxi,maxr,maxc,sub1,sub2]=...
%              BLOCKCORR(im1,im2,[cr cc],[er ec],[cor coc],[cer cec])
%
% "Block-correlate" two "images", or rather [mxnx1] matrices, by taking
% one "tile" from the first and calculating the "correlation coefficient"
% over a range of shifts about the second. The result is the best
% possible position of the particular tile in the first image into the
% second, based on the correlation coefficient. 
%
% INPUT
%
% im1        The first image [mxnx1]
% im2        The second image, of exactly the same size
% cr,cc      Row and column index of center of tile in im1
% er,ec      Row and column extent of tile in first image
% cor,coc    Row and column index of center of correlation search
% cer,cec    Row and column extent of correlation search (could be 0)
%
% The subimage always has odd dimensions, so the specified center is
% really at the center. If you're expecting slight shifts, make cr and cc
% equal to cor and coc and make cer and cec small relative to er and ec.
%
% OUTPUT
%
% sr,sc     Row/col shift (plus is down and right) to im2 to align with im1
% corm      The value of the correlation coefficient at the best shift
% corv      The matrix with the correlation coefficients at shifts
% maxi      Index into corv where we have located the maximum
% maxr,maxc Row and column of im2 to achieve best fit with im1 
% sub1      The "template" out of the im1 (being correlated)
% sub2      The best fitting tile out of im2 (being evaluated)
%
% EXAMPLE  
%
% blockcorr('demo1')
% blockcorr('demo2')
% blockcorr('demo3')
%
% Last modified by fjsimons-at-alum.mit.edu, April 30th, 2008

if ~isstr(im1) % Not the demo
  % Supply defaults
  defval('crcc',[500 600])
  defval('erec',[256 256])
  defval('corcoc',[500 600])
  defval('cercec',[20 30])

  % Extract inputs
  [cr,cc]=deal(crcc(1),crcc(2));
  [er,ec]=deal(erec(1),erec(2));
  [cor,coc]=deal(corcoc(1),corcoc(2));
  [cer,cec]=deal(cercec(1),cercec(2));

  % Extract the relevant tile out of the first image
  hr=floor(er/2);
  hc=floor(ec/2);
  sub1=im1(cr-hr:cr+hr,cc-hc:cc+hc);

  % Define the correlation center positions in the second image
  hcr=floor(cer/2);
  hcc=floor(cec/2);
  [ccr,ccc]=meshgrid(cor-hcr:cor+hcr,coc-hcc:coc+hcc);
  ccr=ccr(:);
  ccc=ccc(:);

  % Initialize the correlation vector
  corv=zeros(length(ccr),1);

  warning off
  % Perform the correlations for each of these guys
  for index=1:length(ccr)
    % Isolate the relevant portion in the second image
    sub2=im2(ccr(index)-hr:ccr(index)+hr,ccc(index)-hc:ccc(index)+hc);
    % Calculate the correlation
    corv(index)=indeks(corrcoef(sub1(:),sub2(:)),2);
  end
  warning on

  % And the later, find the maximum
  [corm,maxi]=max(corv);
  % The column
  maxr=ccr(maxi);
  maxc=ccc(maxi);
  sub2=im2(ccr(maxi)-hr:ccr(maxi)+hr,ccc(maxi)-hc:ccc(maxi)+hc);

  % Return a properly formatted correlation coefficient
  % Do this for ccr and ccc and see that it works
  corv=reshape(corv,cec+~mod(cec,2),cer+~mod(cer,2))';

  % Return the row and column shift to im2 needed to align with im1
  sr=cr-maxr;
  sc=cc-maxc;

  % Output
  vars={sr,sc,corm,corv,maxi,maxr,maxc,sub1,sub2};
  varargout=vars(1:nargout);

elseif strmatch(im1,'demo1')
  [im1,im2]=getimg;
  % Do a first test 
  [sr,sc,cm,cv,mi,mr,mc,s1,s2]=blockcorr(im1,im2,...
			     [500 600],[256 256],[500 600],[20 30]);
  [ah,ha]=krijetem(subnum(2,3));
  axes(ha(1))
  imagesc(im1); axis image
  axes(ha(2))
  imagesc(im2); axis image
  axes(ha(3))
  imagesc(s1); axis image
  axes(ha(4))
  imagesc(s2); axis image
  axes(ha(5))
  imagesc(cv); axis image
  title(sprintf('Center coordinates of image 1 in image 2: %i, %i',...
		mr,mc))

  % Do a first test which should come up with the same answer
  [sr,sc,cm,cv,mi,mr,mc,s1,s2]=blockcorr(im1,im2,...
			     [500 600],[256 256],[500 600],[34 36]);
  axes(ha(6))
  imagesc(cv); axis image
  title(sprintf('Center coordinates of image 1 in image 2: %i, %i',...
		mr,mc))
elseif strmatch(im1,'demo2')
  [im1,im2]=getimg;
  % Do a couple more tests which should all come up with the same answer
  [sr,sc,cm,cv,mi,mr,mc,s1,s2]=blockcorr(im1,im2,...
			     [125 180],[200 200],[125 180],[40 40]);
  [ah,ha]=krijetem(subnum(2,2));
  axes(ha(1))
  imagesc(s1); axis image; title('Image tile 1')
  axes(ha(2))
  imagesc(s2); axis image; title('Image tile 2')
  axes(ha(3))
  imagesc(s2); axis image; title('Image tile 2')
  axes(ha(4))
  imagesc(s1); axis image; title('Image tile 1')
  colormap(gray(256))
  display(sprintf('Answer   is      %i and %i',sr,sc))
elseif strmatch(im1,'demo3')
  [im1,im2]=getimg2;
  % Do a couple more tests which should all come up with the same answer
  [sr,sc,cm,cv,mi,mr,mc,s1,s2]=blockcorr(im1,im2,...
			     [125 180],[200 200],[125 180],[40 40]);
  [ah,ha]=krijetem(subnum(2,2));
  axes(ha(1))
  imagesc(s1); axis image; title('Image tile 1')
  axes(ha(2))
  imagesc(s2); axis image; title('Image tile 2')
  axes(ha(3))
  imagesc(s2); axis image; title('Image tile 2')
  axes(ha(4))
  imagesc(s1); axis image; title('Image tile 1')
  colormap(flipud(gray(2)))
  display(sprintf('Answer   is      %i and %i',sr,sc))
end

% Demo function - get two real images
function [im1,im2]=getimg 
% Get two images
diro='/home/fjsimons/STUDENTS/CatherineRose';
a=ls2cell(fullfile(diro,'*bmp'));
redit=imread(fullfile(diro,a{1}));
[m,n,o]=size(redit);
reds=nan(m,n,2);
greens=nan(m,n,2);
blues=nan(m,n,2);
greys=nan(m,n,2);

for in=1:length(a)
  if in>1
    redit=imread(fullfile(diro,a{in}));
  end
  reds(:,:,in)=reshape(indeks(redit,1:m*n),m,n);
  greens(:,:,in)=reshape(indeks(redit,[1:m*n]+m*n),m,n);
  blues(:,:,in)=reshape(indeks(redit,[1:m*n]+2*m*n),m,n);
  greys(:,:,in)=(reds(:,:,in)+blues(:,:,in)+greens(:,:,in))/3;
end
im1=greys(:,:,1);
im2=greys(:,:,2);

% Demo function - make two fake images
function [im1,im2]=getimg2
% Make two images
im1=ones(1985,1778);
im2=ones(1985,1778);
i1=100; j1=100;
im1(i1,j1)=10;
i2=112; j2=93;
im2(i1,j2)=10;
display(sprintf('Answer should be %i and %i',i1-i2,j1-j2))

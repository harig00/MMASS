function  [mat,c11m,cmnm]=blockisolate2(matrix,c11,cmn,xp,yp,bx,by)
% [mat,c11m,cmnm]=BLOCKISOLATE2(matrix,c11,cmn,xp,yp,bx,by)
%
% Isolates a block of size [bx by] from a matrix with physical
% coordinates c11 and cmn at the center location [xp yp]
%
% INPUT:
%
% matrix         Data field, an MxN matrix
% c11            Physical coordinates of matrix element (1,1)
% cmn            Physical coordinates of matrix element (M,N)
% xp             X-coordinates of the center of the requested block tile
% yp             Y-coordinates of the center of the requested block tile
% bx             Tile dimension, number of samples in the X-direction
% by             Tile dimension, number of samples in the Y-direction
%
% OUTPUT:
%
% mat            The requested block
% c11m           The top left physical coordinates of this block
% cmnm           The bottom right physical coordinates of this block
%
% EXAMPLE:
%
% blockisolate2('demo');
%
% SEE ALSO:
%
% BLOCKISOLATE
%
% Last modified by fjsimons-at-alum.mit.edu, 04/30/2008

if ~isstr(matrix)

  m=size(matrix,1);
  n=size(matrix,2);

  [ind,colnr,rownr]=cor2ind(xp,yp,c11,cmn,m,n);

  mat=matrix(rownr-bx/2:rownr+bx/2-1,colnr-by/2:colnr+by/2-1);

  indo=sub2ind([m n],[rownr-bx/2 rownr+bx/2-1],[colnr-by/2 colnr+by/2]);

  [lo,la]=ind2cor(indo,c11,cmn,m,n);

  c11m=[lo(1) la(1)];
  cmnm=[lo(2) la(2)];
elseif strmatch(matrix,'demo')
  ddir= '/home/fjsimons/MyPapers/2003/JGR-2003/DATA/';
  load(fullfile(ddir,'tdint'))
  clf
  subplot(211)
  C11=[112.6989 -5.7875];
  CMN=[156.9554 -48.2125];
  imagef(C11,CMN,tdint); axis image
  cb=cax2dem([-7400 1500]);
  load(fullfile(getenv('IFILES'),'STATIONS','stations'))
  stt=struct2cell(stations);
  stt=[stt{:}]; stt=reshape(stt,2,length(stt)/2)';
  % Select only Australian stations - ad hoc fix
  stt=stt(stt(:,1)>C11(1) & stt(:,1)<CMN(1)...
	  & stt(:,2)>CMN(2) &  stt(:,2)<C11(2),:)
  [flo,tla]=project(stt(:,1),stt(:,2),C11,CMN);
  hold on; sjuf=shuffle([flo tla]);
  nmsh=5;
  ps=plot(sjuf(1:nmsh,1),sjuf(1:nmsh,2),'v');
  flo=sjuf(1:nmsh,1); tla=sjuf(1:nmsh,2);
  [on,tw,XY]=plotcont([90 10],[180 -60]); delete(tw)
  [clon,clat]=project(XY(:,1),XY(:,2),C11,CMN);
  subplot(212)
  plot(clon,clat,'r'); axis image; axis([C11(1) CMN(1) CMN(2) C11(2)])
  hold on; 
  for index=1:nmsh
    [mat(:,:,index),c11m(index,:),cmnm(index,:)]=...
	blockisolate2(tdint,C11,CMN,flo(index),tla(index),128,128);
    imagef(c11m(index,:),cmnm(index,:),mat(:,:,index))
  end
  ps2=plot(flo(1:nmsh),tla(1:nmsh),'v');
  set([ps ps2],'markerf','w','markere','k')
  cb=cax2dem([-7300 1500]);
end


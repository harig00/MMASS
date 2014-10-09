function varargout=deltaJl(QUAKES,L,method,...
			   rad,nn,el,ww,...
			   U,V,P,dUdr,dVdr,...
			   rho);
% [deltaJ,QUAKES]=deltaJl(QUAKES,L,method,rad,nn,el,ww,U,V,dUdr,dVdr,rho)
%
% Computes the temporal changes in the geodetic coefficients J_l
% due to the earthquake moment tensor in QUAKES.
% 
% INPUT:
%
% QUAKES         Output of READCMT containing earthquake parameters:
%                [time depth lat lon Mtensor]
% L              Spherical harmonic degrees at which J_l is to be computed
% method         1 Explicit method suitable for low-degree computation
%                2 Computation via EQPOTENTIAL, which should do the same thing 
% rad,nn,el,ww   Mode identifying output from GETSPHEROIDAL
% U,V,dUdr,dVdr  Mode eigenfunctions also from GETSPHEROIDAL
% rho            Density profile output from GETMODEL
%
% OUTPUT:
%
% deltaJ         A matrix containing the changes to Jl. Each
%                is a particular event. The first row contains
%                the serial date number of that event, each of
%                the following rows correspond to the order of L.
%                No matter how L is input, the output, beginning
%                with the second row, is in ascending order 
% QUAKES         Output, should you want it
%
%
% EXAMPLES:
%
% deltaJl('demo1',1) % Attempts to recreate Chao & Gross 1987 figure
% deltaJl('demo1',2) % Same thing, with method 2
% deltaJl('demo3') % Same thing, with method 1, until today
% deltaJ=deltaJl(readCMT('demo3'))
%
% Last modified by efwelch-at-princeton.edu, 06/28/2010
% Last modified by fjsimons-at-alum.mit.edu, 07/28/2010

% Check MATLABPOOL and PARFOR to parallelize these loops here

defval('QUAKES',[])
defval('xver',0)

if ~isstr(QUAKES)
  defval('method',1)
  defval('L',[2 3 4 5])
  defval('QUAKES',readCMT('demo1'));
    
  if nargin<4
    % Get the spheroidal mode eigenfunctions
    [rad,nn,el,ww,U,V,P,dUdr,dVdr]=getspheroidal;
    
    % Restrict this to the degrees needed
    ix=find(el>=min(L) & el<=max(L));
    el=el(ix); nn=nn(ix); ww=ww(ix); 
    U=U(:,ix); V=V(:,ix); P=P(:,ix);
    dUdr=dUdr(:,ix); dVdr=dVdr(:,ix);
		    
    % Get the radii and densities of the Earth model
    [radmod,rho]=getmodel;

    % Radius of Earth [m]
    radius=max(rad);
    nrad=length(rad);
    
    % Check these radii are close enough
    difer([rad-radmod]/radius,5,[],NaN)
  end

  % Mass of Earth
  Mass=4*pi*trapeze(rad,rho.*rad.^2);
  
  % Initialize the output
  nQUAKES=size(QUAKES,1);
  deltaJ=nan(1+length(L),nQUAKES);

  % Store times in the deltaJ array
  deltaJ(1,:)=QUAKES(:,1);
  
  % Ensure L is an ordered vector in the return
  L=sort(L(:));
  
  % Restrict to m=0 because this is for Jl only
  % Though in principle the below should work for m>0 also
  m=0;
  
  switch method
   case 1 % Directly computing the strain tensor
    % Get logicals that indicate where each degree is stored
    for j=1:length(L) 
      elis(:,j)=el==L(j);
    end
    % How about the alternative
    % elis=~[repmat(el,1,length(L))-repmat(L(:)',size(el,1),1)];

    % How many unique branch numbers do we have?
    nlen=length(unique(nn));
    
    % Initialize temporary eigenfunction matrices
     Utmp  =zeros(nrad,nlen,length(L));
     Vtmp  =zeros(nrad,nlen,length(L));
    dUdrtmp=zeros(nrad,nlen,length(L));  
    dVdrtmp=zeros(nrad,nlen,length(L));
      wwtmp=zeros(nlen,length(L));
      
    % Assign values for the needed eigenfunctions
    for j=1:length(L);
      Utmp(:,1:sum(elis(:,j)),j)=U(:,elis(:,j));  
      Vtmp(:,1:sum(elis(:,j)),j)=V(:,elis(:,j));  
      dUdrtmp(:,1:sum(elis(:,j)),j)=dUdr(:,elis(:,j));      
      dVdrtmp(:,1:sum(elis(:,j)),j)=dVdr(:,elis(:,j));  
      wwtmp(1:sum(elis(:,j)),j)=ww(elis(:,j));
    end
    
    % Initialize F_nl 
    Fnl=zeros(length(L),nlen);
    % Compute F_nl (C&G eq. 29)
    for j=1:length(L)
      Fnl(j,:)=L(j)*wwtmp(:,j).^(-2).*...
	     trapeze(rad,repmat(rho.*rad.^(L(j)+1),1,nlen).*...
	       (Utmp(:,:,j)+(L(j)+1).*Vtmp(:,:,j)));
    end
    
    % Remove any NaNs that come from dividing by zero
    Fnl(isnan(Fnl))=0;
    
    % Loop over the earthquakes
    h=waitbar(0,sprintf('Looping over all %i earthquakes, method %i',...
			nQUAKES,method));
    for i=1:nQUAKES
      % Get the radius of the hypocenters, at the source [m]
      rs=radius-QUAKES(i,2).*1000;
      
      % Get colatitude (convert to radians), radius of source
      theta=(90-QUAKES(i,3))*pi/180;
      
      % Convert moment tensor [dyne cm] to Joules [J]
      QUAKES(i,5:10)=QUAKES(i,5:10)*1e-7;

      % Find (index of) rad nearest rs at the source
      rsint=findrad(rad,rs,3); rs=rad(rsint);
      
      % Restrict eigenfunctions and derivatives to rs
       Urs  =squeeze(Utmp(rsint,:,:));
      dUdrrs=squeeze(dUdrtmp(rsint,:,:));
       Vrs  =squeeze(Vtmp(rsint,:,:));
      dVdrrs=squeeze(dVdrtmp(rsint,:,:));

      % Later on we expect this to behave like this:
      if prod(size(Urs))==length(Urs)
   	 Urs  =Urs(:);     Vrs  =Vrs(:);
	dUdrrs=dUdrrs(:); dVdrrs=dVdrrs(:);
      end
      
      % Compute spherical harmonics and their theta derivatives
      % This seems to work for positive nonzero orders also
      if L(1)==2 && all(diff(L)==1) && max(L)<15
	% The cumprod only works if L is sequential from 2, if not you need
	% to do factorials... or rather the option after "else" below
	zP1=plm(L,m,cos(theta)); zP2=plm(L-1,m,cos(theta));
	Ylm =(-1)^m*sqrt((2*L+1).*cumprod(L-m)./(4*pi*(cumprod(L+m)))).*zP1;
	Ylmp=(-1)^m*sqrt((2*L+1).*cumprod(L-m)./(4*pi*(cumprod(L+m))))...
	     .*(L*cot(theta).*zP1-(L+m)*csc(theta).*zP2);
      else
	for ind=1:length(L)
	  [Ylem(ind,1),Ydotlem(ind,1)]=...
	      libbrecht(L(ind),cos(theta),'fnc',[],m);
	  Ylm(ind)=Ylem(ind)*sqrt(2-[m==0]);
	  Ylmp(ind)=Ydotlem(ind)*sqrt(2-[m==0]);
	end
      end

      if xver==1
       % Check using computations with YLM, LIBBRECHT, LEGENDREDIFF
       % This seems to work for positive nonzero orders also
       for ind=1:length(L)
	 difer(Ylm(L==L(ind))-ylm(L(ind),m,theta,pi/2),[],[],NaN);
	 [Ylem,Ydotlem]=libbrecht(L(ind),cos(theta),'fnc',[],m);
	 difer(Ylm (L==L(ind))-sqrt(2-[m==0])*Ylem   ,[],[],NaN)
	 difer(Ylmp(L==L(ind))-sqrt(2-[m==0])*Ydotlem,[],[],NaN)
	 if m==0
	   [Ydotl0,~,Yl0]=legendrediff(L(ind),cos(theta),'fnr');
	   difer(Ylm (L==L(ind))-Yl0   /sqrt(4*pi),[],[],NaN)
	   difer(Ylmp(L==L(ind))+Ydotl0/sqrt(4*pi)*sin(theta),[],[],NaN)
	 end
       end
       disp(sprintf('All %i cases verified',ind))
      end
      
      % Compute strain elements for Jl (which means spherical harmonic l=L, m=0)
      % If size(L)==1 then squeeze took out a different singleton
      % dimension in which case you end up having trouble
      E1=dUdrrs.*repmat(Ylm',nlen,1);
      E2=(Urs.*repmat(Ylm',nlen,1)-Vrs.*...
	  repmat((Ylmp.*cot(theta)+L.*(L+1).*Ylm)',nlen,1))/rs;
      E3=(Urs.*repmat(Ylm',nlen,1)+Vrs.*...
	 repmat((Ylmp.*cot(theta))',nlen,1))/rs;
      E4=repmat(Ylmp',nlen,1).*(dVdrrs+(Urs-Vrs)/rs);
      % Obviously this will be zero for m=0, and watch the complexity!
      E5=m*sqrt(-1)*repmat(Ylm',nlen,1)*csc(theta).*...
	 (dVdrrs+(Urs-Vrs)/rs);
      E6=m*2*sqrt(-1)*Vrs*csc(theta)/rs.*...
	 repmat((Ylmp-Ylm*cot(theta))',nlen,1);
      % Perform the mode summation
      vectorsum=[diag(Fnl*E1) diag(Fnl*E2) diag(Fnl*E3) diag(Fnl*E4)...
		 diag(Fnl*E5) diag(Fnl*E6)];
      % Return the deltaJl by the multiplication with the moment tensor
      % If nonzero order, the real/imaginary parts will be Clm and Slm
      deltaJ(2:end,i)=vectorsum*QUAKES(i,5:10)';
      
      % Report on progress
      waitbar(i/nQUAKES)
    end
    % Incoporate scale factor
    deltaJ(2:end,:)=deltaJ(2:end,:).*...
	repmat(-sqrt(2*L+1)*2*sqrt(pi)./...
	       ((2*L+1).*radius.^L*Mass),1,nQUAKES);
   case 2 % By calling EQPOTENTIAL
    % Loop over the earthquakes
    h=waitbar(0,sprintf('Looping over all %i earthquakes, method %i',...
			nQUAKES,method));
    for i=1:nQUAKES
      % Get the radius of the hypocenters, at the source [m]
      rs=radius-QUAKES(i,2)*1000;

      % Compute the strain eigenfunctions nearest the source
      Sstrain=smodestrain(rs,[],[],rad,nn,el,ww,U,V,P,dUdr,dVdr);

      % Calculate the spherical harmonic expansion of the potential
      % Note that this is without further scaling in units of potential
      % that can be plotted and expanded using PLM2XYZ.
      [phiE1,phiE1r]=eqpotential(QUAKES(i,5:10),QUAKES(i,[2 4 3]),...
				 max(L),Sstrain);
      
      % If we are only looking at the m=0 terms, ANY longitude will do!
      % And clearly this works for any positive m
      phiE1r=phiE1r(phiE1r(:,2)==m,:);
      
      % Assign to the correct spot 
      for j=1:length(L);
	% For Clm take column 3, for Slm column 4
	deltaJ(1+j,i)=phiE1r(phiE1r(:,1)==L(j),3);
      end
      
      % Report on progress
      waitbar(i/nQUAKES)
    end

    % Incoporate scale factor that returns the Chao & Gross convention
    % which pulls out the GM/r factor, and note that my PLM2XYZ and
    % PLM2POT are for 4*pi normalized harmonics (for all m and regardless
    % of it). And now I make my coefficients bigger so that THEY can
    % expand in smaller unit-normalized harmonics, as Chao & Gross. 
    deltaJ(2:end,:)=deltaJ(2:end,:)/fralmanac('GravCst')/Mass*radius...
	*sqrt(4*pi);
  end

  close(h)

elseif strcmp(QUAKES,'demo1')
  defval('L',[])
  defval('method',L)
  defval('method',1)
  % do the calculation
  L=[2 3 4 5];

  % Load it or make it and save it
  fname=fullfile(getenv('IFILES'),'CMT',...
		 sprintf('deltaJl_demo1_%i.mat',method));
  if exist(fname,'file')==2 
    load(fname)
    disp(sprintf('Loading %s',fname))
  else
    [deltaJ,QUAKES]=deltaJl([],L,method);
    save(fname,'deltaJ','QUAKES')
  end

  % x-axis limits
  yr1='1976/07/30';
  yr2='1985/06/30';

  ylims=minmax(cumsum(deltaJ(2:end,:),2));
  
  % Just be very explicit to compare the two methods
  ylims=[-0.182964639558030   0.110282156699622]*1e-11;  
  ylims=[-0.224153862361077   0.129030977375356]*1e-11;

  % Ticks
  yrt{1}='1977/01/01'; yrt{2}='1979/01/01'; yrt{3}='1981/01/01';
  yrt{4}='1983/01/01'; yrt{5}='1985/01/01';
  yrtl={'1977','1979','1981','1983','1985'};

  % Make and print the plot
  makeplots(deltaJ,L,yr1,yr2,ylims,yrt,yrtl)
  figdisp([],sprintf('Chao+87_%i',method),[],1) 
elseif strcmp(QUAKES,'demo3')
  defval('L',[])
  defval('method',L)
  defval('method',1)
  % do the calculation
  L=[2 3 4 5];

  % Load it or make it and save it
  fname=fullfile(getenv('IFILES'),'CMT',...
		 sprintf('deltaJl_demo3_%i.mat',method));
  if exist(fname,'file')==2 & 1==3
    load(fname)
    disp(sprintf('Loading %s',fname))
  else
    [deltaJ,QUAKES]=deltaJl(readCMT('demo3'),L,method);
    save(fname,'deltaJ','QUAKES')
  end

  % x-axis limits
  yr1='1976/07/30';
  yr2='2010/09/28';

  ylims=minmax(cumsum(deltaJ(2:end,:),2));

  % Make and print the plot
  makeplots(deltaJ,L,yr1,yr2)
  figdisp([],sprintf('Chao+87_%i',method),[],1) 
end     

% Prepare output
varns={deltaJ,QUAKES};
varargout=varns(1:nargout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeplots(deltaJ,L,yr1,yr2,ylims,yrt,yrtl)

defval('ylims',[])
defval('yrt',[])
defval('yrtl',[])

% Make some plots
clf
[ah,ha,H]=krijetem(subnum(2,2));
cumtime=deltaJ(1,:);

% Single or double axes?
dblax=1;
for ind=1:length(L)
  toplot=cumsum(deltaJ(ind+1,:));

  axes(ah(ind))
  p(ind)=plot(cumtime,toplot);
  if ~isempty(ylims)
    % All axes same limits if you input it
    set(ah(ind),'ylim',ylims*1.075)
  end
  % And all axes same tick marks as the first one
  ytix=get(ah(1),'ytick');
  if dblax==1
    delete(p(ind))
    fax=sqrt(2*L(ind)+1);
    [ax(ind,:),p(ind),p2(ind)]=...
	plotyy(cumtime,toplot,cumtime,toplot/fax);
    if ~isempty(ylims)
      set(ah(ind),'ylim',ylims*1.075)
      set(ah(ind),'ytick',ytix)
    end
    % set(ax(ind,2),'ylim',ylims*1.075/fax)
    set(ax(ind,2),'ylim',get(ax(ind,1),'ylim')/fax)
    % set(ax(ind,2),'ytick',ytix/fax)
    set(ax(ind,2),'ytick',get(ax(ind,1),'ytick')/fax)
    % Rounded yticklabels - figure out the exponent
    ynum=get(ax(ind,2),'ytick');
    ynums=str2num(get(ax(ind,2),'ytickl'));
    warning off MATLAB:divideByZero
    bla=round(log10(ynum(:)./ynums(:)));
    % Verify later that this is well done!
    expo=nanmean(bla(~isinf(bla)));
    disp(sprintf('Exponent on second axis of plot %i is %i',ind,expo))
    warning on MATLAB:divideByZero
    roundyl=num2str(round(ynums*10)/10);
    mspace=str2mat(repmat(32,size(roundyl,1),1));
    % Note that this changes the "auto-exponent" and 'YTickLabelMode'
    set(ax(ind,2),'ytickl',[mspace roundyl mspace])
    axes(ax(ind,2))
    ylx(ind)=ylabel(sprintf('%sC_{l0} (x 1e%i)','\Delta',expo));
    axes(ax(ind,1))
  end
  datetick('x');
  grid on
  xl(ind)=xlabel('year');
  yl(ind)=ylabel(sprintf('%sJ_%d','\Delta',L(ind)));
end     

% Cosmetics
set(ha(3:4),'yaxisl','r')
longticks(ah(:))
set(ah,'xlim',[datenum(yr1) datenum(yr2)])
% Despite the above DATETICK, want to be more explicit
if ~isempty(yrt)
  set(ah,'xtick',datenum(yrt),...
	 'xtickl',yrtl)
end
set([xl yl],'color','k')
delete(xl(1:2))
nolabels(ah(1:2),1)
[bh,th]=label(ah,'ll',[],[],[],[],[],1.05);
serre(H',1/2.25,'down')

if dblax==1
  noticks(ax(:,2),1)
  delete(ylx([2 4]))
  delete(p2)
  set(ax([2 4],2),'yaxisl','l')
  set(ax,'YColor','k')
  set([ax(:,2)' ylx([1 3])],'FontS',8)
  % Was /50
  moveh(ylx([1 3]),range(get(ah(1),'xlim'))/75)
end

% Refer to WGS84
GM=fralmanac('GM_wgs84');
rf=fralmanac('rf_wgs84');
a=fralmanac('a_wgs84');
omega=fralmanac('omega_wgs84');

% Modify the coefficients to the reference ellipsoid
halfL=2;
Jevens=-kindeks(grs(GM,rf,a,omega,halfL),1).*sqrt(2*2*[1:halfL]'+1);

% Add the even-degree Jl factors from WGS-84
for ind=[1 3]
  [bh(ind),th(ind)]=...
      boxtex('ur',ah(ind),sprintf('J_%i = %8.3e',L(ind),...
				  Jevens(round(ind/2))),...
	     10,[],[],1.10,1.25);
  movev(th(ind),-range(get(ah(ind),'ylim'))/50)
end

fig2print(gcf,'portrait')

function varargout=admitos(thhats,th0,params,ah)
% [ah,ha,yl]=ADMITOS(thhats,th0,params,ah)
%
% Plots the admittances as they have been recorded in the result of a
% series of simulations. All of these are consistent with the data being
% given by the "true" admittance and thus amount to an error estimate on
% the admittance given these data. As seen in Olhede & Simons. 
%
% INPUT:
%
% thhats     The estimated model parameter vector
% th0        The true model parameter vector
%            th0(1)=D    Isotropic flexural rigidity 
%            th0(2)=f2   The sub-surface to surface initial loading ratio 
%            th0(3)=s2   The first Matern parameter, aka sigma^2 
%            th0(4)=nu   The second Matern parameter 
%            th0(5)=rho  The third Matern parameter 
% params     A parameter vector with constants assumed known, see SIMULOS
%            In order of apppearance and as needed only
%            DEL   surface and subsurface density contrast [kg/m^3]
%            g     gravitational acceleration [m/s^2]
%            z2    the positive depth to the second interface [m]
% ah         Two axis handles if you already have them
%
% OUTPUT:
%
% ah,ha      Various axis handles of the plots made
%
% EXAMPLE:
%
% This only gets used in MLEOS thus far
%
% SEE ALSO:
%
% MLEOS, MLEPLOS, ADMITTANCE, FORSYTH
%
% Last modified by fjsimons-at-alum.mit.edu, 04/17/2012

% Default values
defval('thhats',[])
defval('th0',[])
defval('params',[])
defval('stdmult',[-2 2])
defval('thth',[10 100 1000]);

% Extract the needed auxiliary parameters in meaningful terms
% The density contrasts [kg/m^3]
DEL=params(1:2);
% The gravitational acceleration [m/s2]
g  =params(3);
% The interface depth [m]
z2 =params(4);

% Extract the true values in meaningful terms
tD =th0(:,1);
tf2=th0(:,2);
ts2=th0(:,3);
tnu=th0(:,4);
trh=th0(:,5);

% Extract the estimated values in meaningful terms
% If you were to plot the histograms you'd get what MLEOS('demo2') does 
D =thhats(:,1);
f2=thhats(:,2);
s2=thhats(:,3);
nu=thhats(:,4);
rh=thhats(:,5);

% Calculate the theoretical admittance/coherence of the truth
% Maybe see if where we evaluate would have been possible with the data
% size? This could be interesting
lbd=logspace(4,6.7,100);
% Young's modulus
defval('E',1.4e11);
% Poisson's ratio
defval('v',0.25);
tTe=DtoTe(tD,E,v);
[tG2b,k,l,tZb,Zf]=forsyth(tTe,lbd,tf2,0,DEL(1),DEL(2),z2,g,1);

% Calculate the theoretical admittance/coherence of the estimates
Te=DtoTe(D,E,v);
[Zb,G2b]=deal(nan(length(k),length(Te)));
for ind=1:length(Te)
  [G2b(:,ind),~,l,Zb(:,ind)]=...
      forsyth(Te(ind),lbd,f2(ind),0,DEL(1),DEL(2),z2,g);
end

% For good measure calculate the transition wavelength
[k12,l12,E,v]=transl(sqrt(tf2),tTe/1000,DEL(1),DEL(2));

% Now the plot
clf
defval('ah',NaN)
if isnan(ah)
  [ah,ha]=krijetem(subnum(1,2));
end

axes(ah(1))
% Plot the results of the truth
ptZb=semilogx(l/1e3,tZb,'k','LineW',1); hold on
set(ah(1),'ylim',[-0.21 0.01])
plot(repmat(thth,2,1),ylim,':')
% Plot the results of the estimation
pZb=semilogx(l/1e3,Zb,'Color',grey,'LineW',0.5); 
top(ptZb,ah(1)); xl(1)=xlabel(1); hold off
yl(1)=ylabel('admittance (mgal/m)');

axes(ah(2))
% Plot the results of the truth
ptG2b=semilogx(l/1e3,tG2b,'k','LineW',1); hold on
p12=semilogx(l12,1/2);
set(p12,'Marker','o','MarkerFaceC','w','MarkerEdgeC','k');
set(ah(2),'ylim',[-0.05 1.05])
plot(repmat(thth,2,1),ylim,':')
% Plot the results of the estimation
pG2b=semilogx(l/1e3,G2b,'Color',grey,'LineW',0.5); 
top(ptG2b,ah(2)); xl(2)=xlabel(1); hold off
top(p12,ah(2));
yl(2)=ylabel('coherence');

% Cosmetics and axis
set(ah(1),'ydir','rev')
serre(ah,[],'across')
set(ah,'xdir','rev','xgrid','off','ygrid','on')
set(ah,'xtick',thth,'xtickl',thth)
set(ah,'xlim',minmax(l)/1000)
set(xl,'string','wavelength (km)')
longticks(ah)

% Now how about this - we know the truth, let's resimulate some data and
% form some quick admittance and coherence plots from them, see how far
% we get using that approach
fields={'DEL','g','z2','dydx','NyNx','blurs'};
defstruct('pars',fields,...
	  {params(1:2),params(3),params(4),params(5:6),params(7:8),params(9)});

G=fralmanac('GravCst');
[k,kx,ky,dci,dcn]=knum2(pars.NyNx,...
		[(pars.NyNx(1)-1)*pars.dydx(1) (pars.NyNx(2)-1)*pars.dydx(2)]);
% Wavenumber in radians per m
kr=k(dci(1),dci(2):end);
% Wavelength in km
lams=2*pi./kr/1e3;

defval('fftplot',1)
defval('fftplots',0)
defval('mtmplot',0)
defval('mtmplots',0)

[coh,adm]=deal(nan(pars.NyNx(2)/2,size(thhats,1)));
for ind=1:length(thhats)*2
  % Since we simulate blurred data there's no need for anything else!
  % We simulate data which when analyzed give the correct result, how
  % good is that! Am I getting this right?
  % When I DO blur I see the edge effects, when I don't I don't  but
  % that's how I want it - we need to "observe" the correctly blurred
  % version of the truth - and yes, we need to acknowledge this blurring
  % later on, but this so we can't get misled. Simply use the defaults.
  clc
  % Blurred or not? If not, you won't see the bias, thinking this is
  % because the windows undo each other, and they shouldn't.

  % BUT HERE YOU SEE HOW SENSITIVE THIS IS TO THE CHOICE OF RHO - WHETHER
  % YOU SEE THE EFFECT OR NOT
  [T,G,th0dd,p,k,Hk,Gk]=simulos(th0,pars,pars.blurs);
  % Coherence and admittance - note position of averaging!
  [~,Ctops]=radavg(v2s([Hk(:,1).*conj(Gk)]));
  [~,Cbot1]=radavg(v2s(Hk(:,1).*conj(Hk(:,1))));
  [~,Cbot2]=radavg(v2s(Gk.*conj(Gk)));
  % Coherence
  coh(:,ind)=abs(Ctops).^2./Cbot1./Cbot2;
  % Admittance in mgal/m
  adm(:,ind)=realize(Ctops./Cbot1/1e-5);
  
  % I'm wondering if Sofia referred to this as needing the mean and
  % variance of lots of numerators/denominators to get an idea of the
  % mean and variance of the coherence as measured by real data. 
  % Rather than the mean and variance of lots of coherences, that
  % wouldn't be right? Or would it? Wouldn't we want to know the
  % distribution of many such estimates as they are actually made? I
  % believe she referred to this type of averaging to come up with a
  % coherence in the first place, though I now do this using azimuthal
  % averaging. I now also need to compare the wavenumber resolution with
  % the NW and the Nyquist!
  
  if mtmplots==1
    NW=3;
    [~,~,~,COH2,~,ADM]=mtm(v2s(T(:,1)),v2s(G),NW);
    [~,COH2]=radavg(COH2);
    [~,ADM]=radavg(ADM/1e-5);
    axes(ah(1)); hold on
    semilogx(lams',ADM,'r+');
    axes(ah(2)); hold on
    semilogx(lams',COH2,'r+');
  end
  
  % Plot each individual estimate
  if fftplots==1
    axes(ah(1)); hold on
    semilogx(lams,adm(:,ind),'r+');
    axes(ah(2)); hold on
    semilogx(lams,coh(:,ind),'r+');
  end
end

defval('conftype',2)

% Plot a summary of the estimations
if fftplot==1
  lamlams=repmat(lams(:),1,2);
  % ADMITTANCE ESTIMATION
  axes(ah(1)); hold on
  switch conftype
    case 1
     am=semilogx(lams,mean(adm,2),'ko');
     set(am,'MarkerE','k','MarkerFaceC','k','MarkerS',3)
     as=semilogx(lamlams',[repmat(mean(adm,2),1,2)+...
     		    [stdmult'*realize(std(adm'))]']','k-');
     case 2
      am=semilogx(lams,median(adm,2),'ko');
      set(am,'MarkerE','k','MarkerFaceC','k','MarkerS',3)
      % Maybe plot 95 confidence interval instead
      as=semilogx(lamlams',prctile(adm',[2.5 97.5]),'k-');
  end

  % COHERENCE ESTIMATION
  axes(ah(2)); hold on
  switch conftype
   case 1
    cm=semilogx(lams,mean(coh,2),'ko');
    set(cm,'MarkerE','k','MarkerFaceC','k','MarkerS',3)
    cs=semilogx(lamlams',[repmat(mean(coh,2),1,2)+...
		    [stdmult'*realize(std(coh'))]']','k-');
    case 2
     cm=semilogx(lams,median(coh,2),'ko');
     set(cm,'MarkerE','k','MarkerFaceC','k','MarkerS',3)
     % Maye plot 95 confidence interval instead
     cs=semilogx(lamlams',prctile(coh',[2.5 97.5]),'k-');
  end
end

% THIS ALL MAKES SENSE NOW, BUT HOW TO PUT IT IN THE PAPER?
  
% So we know the true admittance. And the MLE finds the true admittance.
% Now we compare how unwindowed FFT would see the admittance. Which has
% edge effects. And a high variance. So we can also 
% Using the multitaper windowing method we can estimate a blurred admittance.

% I suppose we should compare this to the theoretical admittance which we
% THEN BLUR. But the shape difference between the curves gives us an idea
% of the variance, we are not that concerned with the systematic bias,
% are we. The variance trumps and dominates the bias. See the paper.

if mtmplot==1
  % Brief peek at single multitaper? Could of course look at the estimate
  % of the power spectral densities also... etc. etc. 
  NW=3;
  [SX,SY,SXY,COH2,vCOH2,ADM]=mtm(v2s(T(:,1)),v2s(G),NW);
  [~,COH2]=radavg(COH2);
  % Variance of the average is not the average of the variance, would need
  % to go all the way back to Dong's and Lara's code to remake sense of
  % this. This would mean an extra division by N as usual
  [~,vCOH2,~,~,~,naz]=radavg(vCOH2); 
  % Well, this should only be the case when you use the truly independent
  % number of estimates. Which we did right in Dong's code, but not here.
  %vCOH2=vCOH2/naz;
  % The wavenumber 'resolution' - see ALL_BBB_1D, including the factor of 2
  dk=2*pi*NW/max(pars.NyNx.*pars.dydx)/2;
  vCOH2=vCOH2/range(minmax(kx(kx>0)))*dk;
  % This is more or less acceptable
  % So we will use this here
  [~,ADM]=radavg(ADM);
  % If this variance "works" should get this everywhere inside the loop,
  % which we have verified.
  axes(ah(1)); 
  semilogx(lams,real(ADM)*1e5,'bo-');
  hold off
  axes(ah(2)); 
  bl=semilogx(lams,COH2,'bo-');
  bk=semilogx(2*pi./[kr(:) kr(:)]'/1e3,[COH2 COH2]'+[sqrt(vCOH2)*stdmult]','b-');
  %bg=semilogx(2*pi./[kr(:)-dk kr(:)+dk]'/1e3,[COH2 COH2]','g-','LineW',2);
  hold off

  % Now we see how FFT has lots of variance but low bias. Multitaper has
  % lots of bias but low variance. BUT EVEN THEN we cannot do a proper
  % estimation for the reasons that we have already understood.
  
  % PERIODOGRAM IS NOT CONSISTENT, VARIANCE DOESN'T DECREASE WITH MORE K
end

% Reorder what's on top
bottom(am,ah(1))
bottom(as,ah(1))
bottom(cm,ah(2))
bottom(cs,ah(2))

% Do this so the reduction looks slightly better
set(yl,'FontS',12)

% Output
varns={ah,ha,yl};
varargout=varns(1:nargout);

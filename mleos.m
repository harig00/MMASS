function varargout=mleos(Hx,Gx,thini,params,algo,bounds,aguess)
% [thhat,covh,logli,thini,scl,params,eflag,oput,grd,hes,Hk,k,options,bounds]=...
%          MLEOS(Hx,Gx,thini,params,algo,bounds,aguess)
%
% Performs a maximum-likelihood estimation for UNCORRELATED loads as in
% Olhede & Simons (2013) by minimization using FMINUNC/FMINCON.
% This code uses an early, abandoned, version of the initial-load scaling.
%
% INPUT:
%
% Hx       Real matrix with surface and subsurface topography [m m], 
% Gx       Real vector with Bouguer gravity anomaly data [m/s^2]
%          ... either from direct observation or from SIMULOS...
%          ... with the geographical coordinates linearly unwrapped...
% thini    A starting guess for the parameter vector with elements:
%          [D f2 s2 nu rho], in Nm, and "nothing", see SIMULOS
% params   A parameter structure with constants assumed known, see SIMULOS
%          [DEL g z2 dydx NyNx blurs kiso] in the units of ...
%           kg/m^3 (2x), m/s^2, m (3x), "nothing" (3x), rad/m
%           blurs  0 Don't blur likelihood using the Fejer window
%                  N Blur likelihood using the [default: N=2] resampled Fejer window
%           kiso   wavenumber beyond which we are not considering the likelihood
% algo     'unc' uses FMINUNC
%          'con' uses FMINCON with positivity constraints [default]
%          'klose' simply closes out a run that got stuck
% bounds    A cell array with those positivity constraints [defaulted]
% aguess    A parameter vector [s2 nu rho] that will be used in
%           simulations for demo purposes, and on which "thini" will be
%           based if that was left blank. If "aguess" is blank, there is
%           a default. If "thini" is set, there is no need for "aguess"
%
% OUTPUT:
%
% thhat    The maximum-likelihood estimate of the vector [scaled]:
%          [D f2 s2 nu rho], in Nm, and "nothing", see SIMULOS
% covh     The asymptotic covariance matrix of this estimate
% logli    The maximized value of the likelihood
% thini    The starting guess used in the optimization
% scl      The scaling applied as part of the optimization procedure
% params   The known constants used inside, see above
% eflag    The exit flag of the optimization procedure [bad if 0]
% oput     The output structure of the optimization procedure
% grd      The gradient of the misfit function at the estimate
% hes      The Hessian of the misfit function at the estimate
% Hk       The spectral-domain interface topographies after deconvolution 
% k        The wavenumbers on which the estimate is actually based
% options  The options used by the optimization procedure
% bounds   The bounds used by the optimization procedure
%
% NOTE: 
%
% At least 'demo1' has been tested to run in an SPMD loop! Files are
% opened in append mode - except for thzro, which will only reflect one lab.
%
% EXAMPLE:
%
%% Perform a series of N simulations centered on th0
% mleos('demo1',N,th0)
%% Statistical study of a series of simulations done using 'demo1'
% mleos('demo2','07-Jan-2013-blurs-2')
%% Admittance/coherence study of a series of simulations
% mleos('demo3','07-Jan-2013-blurs-2')
%% Covariance study of a series of simulations
% mleos('demo4','07-Jan-2013-blurs-2')
%% One simulation and a chi-squared plot
% mleos('demo5')
%
% Last modified by fjsimons-at-alum.mit.edu, 02/13/2015

if ~isstr(Hx)
  defval('algo','con')
  if strcmp(algo,'con')
    % Parameters for FMINCON in case that's what's being used
    bounds={[],[],... % Linear inequalities
	    [],[],... % Linear equalities
	    [1e17 eps eps 0.95/100 1]*100,... % Lower bounds
	    [Inf Inf Inf 4.00     Inf],... % Upper bounds
	    []}; % Nonlinear (in)equalities
  else
    bounds=[];
  end
  % The necessary strings for formatting
  str0='%27s';
  str1='%12.0e ';
  str2='%12.5g ';
  
  % Supply the needed parameters, keep the givens, extract to variables
  fields={'DEL','g','z2','dydx','NyNx','blurs','kiso'};
  defstruct('params',fields,...
	    {[2670 630],9.81,35000,[20 20]*1e3,sqrt(length(Hx))*[1 1],2,NaN});
  struct2var(params)
  
  % The gravitational constant (in m^3/kg/s2)
  G=fralmanac('GravCst');

  % Being extra careful or not?
  defval('xver',0)
  
  % The parameters used in the simulation for demos, or upon which to base "thini"
  defval('aguess',[1e24 0.8 0.0025 2 2e4]);
  % Scale the parameters by this factor; fix it unless "thini" is supplied
  defval('scl',10.^round(log10(abs(aguess))));

  % Unless you supply an initial value, construct one from "aguess" by perturbation
  nperturb=0.25;
  defval('thini',abs((1+nperturb*randn(size(aguess))).*aguess))
  disp(sprintf(sprintf('%s : %s ',str0,repmat(str2,size(thini))),...
	       'Starting theta',thini))

  % If you brought in your own initial guess, need an appropriate scale
  if ~isempty(inputname(3)) || any(aguess~=thini)
    scl=10.^round(log10(abs(thini)));
    disp(sprintf(sprintf('%s : %s ',str0,repmat(str1,size(scl))),...
		 'Scaling',scl))
  end
  % Now scale so the minimization doesn't get into trouble - bounds also!
  thini=thini./scl;
    
  defval('taper',0)
  if taper==1
    % Were going to want to make a 2D taper - any taper
    disp(sprintf('%s with TAPERING, DO NOT DO THIS YET',upper(mfilename)))
    NW=2;
    E1=dpss(NyNx(1),NW,1);
    E2=dpss(NyNx(2),NW,1);
    Tx=repmat(E1,1,NyNx(2)).*repmat(E2',NyNx(1),1);
    % But should still watch out for the spectral gain I suppose, this isn't
    % done just yet, although it appears already properly normalized
    % However, this looks better now, doesn't it?
    Tx=Tx*sqrt(prod(NyNx));
    % Not doing anything still amounts to saying Tx=1
  else
    Tx=1;
  end

  % Create the appropriate wavenumber axis
  k=knums(params);

  % Modify to demean
  disp('DEMEAN BOTH DATA SETS')
  Hx(:,1)=Hx(:,1)-mean(Hx(:,1));
  Gx=Gx-mean(Gx);
  % Let us NOT demean and see where we end up...

  % Turn the observation vector to the spectral domain
  Hk(:,1)=tospec(Tx(:).*Hx(:,1),params);
  Gk     =tospec(Tx(:).*Gx     ,params);
  Hk(:,2)=Gk.*exp(k(:).*z2)/2/pi/G/DEL(2);

  if xver==1
    % This is reliant on the gravity data to be amenable to deconvolution,
    % while we could fake it in the simulations by working with Hx. Check
    % quickly, and note that there are roundoff errors right away!
    % Is the normalization right? I recently absorbed this into TOSPEC.
    difer([tospec(Tx(:).*Hx(:,2),p)-Hk(:,2)]/length(Hk),8,[],NaN)
    % Should also compare this with what actually can come out of SIMULOS
    % itself, although with real data of course we don't have this.
  end
  
  NN=200;
  % And now get going with the likelihood using Hk(:,1:2) or [Hk(:,1) Gk]
  % [ off | iter | iter-detailed | notify | notify-detailed | final |
  % final-detailed ] 
  % Should probably make the tolerances relative to the number of k points
  options=optimset('GradObj','off','Display','final',...
		   'TolFun',1e-11,'TolX',1e-11,'MaxIter',NN,...
		   'LargeScale','off');
  % The 'LargeScale' option goes straight for the line search when the
  % gradient is NOT being supplied.

  % Set the parallel option to (never) use it for the actual optimization
  % Doesn't seem to do much when we supply our own gradient
  options.UseParallel='always';

  if blurs==0 || blurs==1
    % Use the analytical gradient in the optimization, rarely a good idea
    % options.GradObj='on';
    if xver==1
      % Definitely good to check this once in a while
      options.DerivativeCheck='on';
    end
  end

  % And find the MLE! Work on scaled parameters
  try
    switch algo
     case 'unc'
      % disp('Using FMINUNC for unconstrained optimization of LOGLIOS')
      t0=clock;
      [thhat,logli,eflag,oput,grd,hes]=...
	  fminunc(@(theta) loglios(theta,params,Hk,k,scl),...
		  thini,options);
      ts=etime(clock,t0);
      % Could here compare to our own estimates of grad and hes!
     case 'con'
      % New for FMINCON
      options.Algorithm='active-set';
      % disp('Using FMINCON for constrained optimization of LOGLIOS')
      t0=clock;
      [thhat,logli,eflag,oput,lmd,grd,hes]=...
	  fmincon(@(theta) loglios(theta,params,Hk,k,scl),...
		  thini,...
		  bounds{1},bounds{2},bounds{3},...
                  bounds{4},bounds{5}./scl,bounds{6}./scl,bounds{7},...
		  options);
      ts=etime(clock,t0);
     case 'klose'
       % Simply a "closing" run to return the options
       varargout=cellnan(nargout,1,1);
       varargout{end-1}=options;
       varargout{end}=bounds;
       return
    end
    if xver==1
      disp(sprintf('%8.3gs per %i iterations or %8.3gs per %i function counts',...
		   ts/oput.iterations*100,100,ts/oput.funcCount*1000,1000))
    else
      disp(sprintf('\n'))
    end
  catch
    % If something went wrong, exit gracefully
    varargout=cellnan(nargout,1,1);
    return
  end
  
  % This is the entire-plane estimate (hence the factor 2!)
  covh=hes2cov(hes,scl,length(k(~~k))*2);

  % Talk!
  disp(sprintf(sprintf('\n%s : %s ',str0,repmat(str2,size(thhat))),...
	       'Estimated theta',thhat.*scl))
  disp(sprintf(sprintf('%s : %s ',str0,repmat(str2,size(thhat))),...
	       'Asymptotic stds',thhat.*scl))

  % Generate output as needed
  varns={thhat,covh,logli,thini,scl,params,eflag,oput,grd,hes,Hk,k,options,bounds};
  varargout=varns(1:nargout);
elseif strcmp(Hx,'demo1')
  % If you run this again on the same date, we'll just add to THINI and
  % THHAT but you will start with a blank THZERO. See 'demo2'
  % How many simulations: the second argument after the demo id
  defval('Gx',[]);
  N=Gx; clear Gx
  defval('N',500)
  % What th-parameter set? The SECOND argument after the demo id
  defval('thini',[]);
  % If there is no preference, then that's OK, it gets taken care of
  th0=thini; clear thini
  % What fixed-parameter set? The THIRD argument after the demo id
  defval('params',[]);

  % The number of parameters to solve for
  np=5;

  % Open files and return format strings
  [fid0,fid1,fid2,fid3,fmt1,fmt2,fmt3,fmtf,fmte,fmtd,fmtc]=...
      osopen(np);
 
  % Do it!
  good=0; 
  % Initialize the average Hessian
  avH=zeros(np,np);

  % Set N to zero to simply close THZERO out
  for index=1:N
    % Simulate data from the same lithosphere, watch the blurring
    [Hx,Gx,th0,p,k]=simulos(th0,params);
    % Check the dimensions of space and spectrum are right
    difer(length(Hx)-length(k(:)),[],[],NaN)

    % Form the maximum-likelihood estimate
    t0=clock;
    [thhat,covh,logli,thini,scl,p,e,o,gr,hs,~,~,ops,bnds]=...
	mleos(Hx,Gx,[],p,[],[],th0);
    ts=etime(clock,t0);

    % Initialize the THZRO file
    if index==1
      oswzerob(fid0,th0,p,ops,bnds,fmt1,fmt2)
    end

    % If a model was found, keep the results, if not, they're all NaNs
    % Ignore the fact that it may be at the maximum number of iterations
    % e=1

    % IF NUMBER OF FUNCTION ITERATIONS IS TOO LOW DEFINITELY BAD
    itmin=0;

    % A measure of first-order optimality (which in this
    % unconstrained case is the infinity norm of the gradient at the
    % solution)  
    % FJS to update what it means to be good - should be in function of
    % the data size as more precision will be needed to navigate things
    % with smaller variance! At any rate, you want this not too low.
    optmin=Inf;
    % Maybe just print it and decide later? No longer e>0 as a condition.
    % e in many times is 0 even though the solution was clearly found, in
    % other words, this condition IS a way of stopping with the solution
    % Remember that the correlation coefficient can be negative or zero!
    % The HS is not always real, might be all the way from the BLUROS?
    try % Because if there are NaNs or not estimate it won't work
	% Maybe I'm too restrictive in throwing these out? Maybe the
	% Hessian can be slightly imaginary and I could still find thhat
	if isreal([logli gr(:)']) ...
	      && all(thhat>0) ...
	      && all(~isnan(thhat)) ...
	      && o.iterations > itmin ...
	      && o.firstorderopt < optmin
	  good=good+1;
	  % Build the average of the Hessians for printout later
	  avH=avH+hs.*[scl(:)*scl(:)'];
	  % Reapply the scaling before writing it out
	  fprintf(fid1,fmt1,thhat.*scl);
	  fprintf(fid2,fmt1,thini.*scl);
	  % Here we compute and write out the moments of the Xk
	  [L,~,momx]=logliosl(thhat,p,Hk,k,scl);
	  
	  % Print the optimization diagnostics to a different file	
	  oswdiag(fid3,fmt1,fmt3,logli,gr,hs,thhat,thini,scl,ts,e,o,....
		  var(Hx),momx,covh)
	end
    end
  end
  % If there was any success at all, finalize the THZERO file
  % If for some reason this didn't end well, do an N==0 run.
  
  % Initialize if all you want is to close the file
  if N==0; 
    [Hx,Gx,th0,p,k]=simulos(th0,params); 
    good=1; avH=avH+1; 
    [~,~,~,~,~,pp,~,~,~,~,~,~,ops,bnds]=mleos(Hx,Gx,[],[],'klose');
    oswzerob(fid0,th0,p,ops,bnds,fmt1,fmt2)
  end
  
  if good>=1
    % This is the scaling based on the truth which we use here 
    sclth0=10.^round(log10(th0));

    % This is the average of the Hessians, should be close to the Fisher
    avH=avH./[sclth0(:)*sclth0(:)']/good;

    % Now compute the theoretical covariance and scaled Fisher
    [covF,F]=covthos(th0./sclth0,p,k,sclth0);
    % Of course when we don't have the truth we'll build the covariance
    % from the single estimate that we have just obtained. This
    % covariance would then be the only thing we'd have to save.
    if labindex==1
      oswzeroe(fid0,sclth0,avH,good,F,covF,fmtc,fmte,fmtf)
    end
  end
  
  % Put both of these also into the thzro file 
  fclose('all');
elseif strcmp(Hx,'demo2')
  defval('Gx',[]);
  datum=Gx;
  defval('datum',date)

  % The number of parameters to solve for
  np=5;

  % Load everything you know about this simulation
  [th0,thhats,params,truecov,E,v,~,~,momx]=osload(datum);

  % Report the findings of the moment parameters
  disp(sprintf('m(m(Xk)) %f m(v(Xk)) %f m(magic) %s v(magic) %f',...
	      mean(momx),var(momx(:,end))))

  % Plot it all - perhaps some selection on optis?
  [ah,ha]=mleplos(thhats,th0,truecov,E,v,params,sprintf('MLEOS-%s',datum));
  
  % Print the figure! Don't forget the degs.pl script
  figna=figdisp([],sprintf('%s_%s',Hx,datum),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna)); 
  system(sprintf('rm -f %s.eps',figna)); 
elseif strcmp(Hx,'demo3')
  defval('Gx',[]);
  datum=Gx;
  defval('datum',date)

  % The number of parameters to solve for
  np=5;

  % Load everything you know about this simulation
  [th0,thhats,params,truecov,E,v]=osload(datum);
  
  % Plot it all: one admittance/coherence curve for every estimate
  % Well, only plot a hundred at random, how about that
  [ah,ha]=admiplos(thhats(randi(length(thhats),100,1),:),th0,truecov,E,v,params,[],length(thhats));
  
  % Make the plot
  figna=figdisp([],sprintf('%s_%s',Hx,datum),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna));
elseif strcmp(Hx,'demo4')
  defval('Gx',[]);
  datum=Gx;
  defval('datum',date)
  
  % The number of parameters to solve for
  np=5;

  % Load everything you know about this simulation
  [th0,thhats,params,truecov,E,v,obscov,sclcov]=osload(datum);
  
  % Make the plot
  ah=covplos(2,sclcov,obscov,truecov,params,thhats,th0,E,v,'ver');

  % Make the plot
  figna=figdisp([],sprintf('%s_%s',Hx,datum),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna)); 
  system(sprintf('rm -f %s.eps',figna)); 
elseif strcmp(Hx,'demo5')  
  % What th-parameter set? The SECOND argument after the demo id
  defval('Gx',[]);
  % If there is no preference, then that's OK, it gets taken care of
  th0=Gx; clear Gx
  % What fixed-parameter set? The THIRD argument after the demo id
  defval('thini',[]);
  params=thini; clear thini
  
  % Figure name
  figna=sprintf('%s_%s_%s',mfilename,Hx,date);

  % Simulate data, watch the blurring, verify COLCHECK inside
  xver=1;
  [Hx,Gx,th0,p,k,Hk]=simulos(th0,params,xver);
  
  % Initialize, take defaulted inside MLEOS for now
  thini=[];

  % Perform the optimization, whatever the quality of the result
  [thhat,~,logli,thini,scl,p,e,o,gr,hs]=mleos(Hx,Gx,thini,p);
  % Maybe here could test using  the "wrong" code, MLEROS or MLEROS0

  % Take a look at the unblurred gradient purely for fun, they should be
  % so small as to be immaterial
  grobs=-nanmean(gammakos(k,thhat.*scl,p,Hk))';
  
  % Take a look at the unblurred theoretical covariance at the estimate,
  % to compare to the observed blurred Hessian; in the other demos we
  % compare how well this works after averaging
  [covthat,F]=covthos(thhat,p,k,scl);

  % Collect the theoretical covariance for the truth for the title
  covth=covthos(th0./scl,p,k,scl);
 
  % Take a look at the scaled Fisher to compare with the scaled Hessian  
  F;
  hs;
  grobs;
 
  % Take a look at the scaled covariances
  predcov=covthat./[scl(:)*scl(:)'];
  % And compare to the inverse of the scaled Hessians in the full plane
  % Watch as this will change with the anticipated halfplane changes
  obscov=inv(hs)/length(k(:))*2;
  % These two should compare reasonably well in the unblurred case, I
  % would have thought - but of course it's ignoring stochastic
  % variability. If we can use the Hessian for the uncertainty estimation
  % we can do this for the cases where we can't come up with a
  % theoretical covariance, not even an unblurred one. Check std's.
  % Maybe should do this a hundred times?
  disp(sprintf('%s\n',repmat('-',1,97)))
  disp('predicted and observed scaled standard deviations and their ratio')
  disp(sprintf([repmat('%6.3f  ',1,length(obscov)) '\n'],...
	       [sqrt(diag(predcov))' ; sqrt(diag(obscov))' ; ...
		sqrt(diag(predcov))'./sqrt(diag(obscov))']'))
  disp(repmat('-',1,97))
  % Talk again!
  [str0,str2]=osdisp(th0,p);
  disp(sprintf(sprintf('%s : %s ',str0,repmat(str2,size(thhat))),...
	       'Estimated theta',thhat.*scl))
  disp(repmat('-',1,97))
  
  % Young's modulus 
  defval('E',1.4e11);
  % Poisson's ratio
  defval('v',0.25);

  % Quick plot, but see OLHEDESIMONS5
  clf
  ah=krijetem(subnum(2,3)); delete(ah(4:6)); ah=ah(1:3);

  % Maybe we should show different covariances than the predicted ones??

  % Time to rerun LOGLIOS one last time at the solution
  [L,~,momx]=loglios(thhat,p,Hk,k,scl);

  % Better feed this to the next code, now it's redone inside
  mlechiplos(1,Hk,thhat,scl,p,ah,0,th0,covth,E,v);

  % Print to file
  figna=figdisp([],sprintf('%s','demo7_2'),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna)); 
  system(sprintf('rm -f %s.eps',figna)); 
end


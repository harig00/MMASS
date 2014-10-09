function varargout=ilk(lemma,l,em,N)
% [errs,l,em,x]=ILK(lemma,l,em,N)
%
% Tests three lemmata by Ilk (1983) quoted by Eshagh (2009a, p. 137),
% modified to be calculated for Schmidt-normalized Legendre associated
% functions as Matlab defines them, without (-1)^m and with the
% sqrt(2-dom), i.e. the way they come out of LEGENDRE and LIBBRECHT 'sch'.
%
% INPUT:
%
% lemma    1 The first  lemma by Ilk (1983, Z.1.44)
%          2 The second lemma by Ilk (1983, Z.1.41)
%          3 The third  lemma by Ilk (1983, Z.1.42)
% l        A particular degree to be tested (random default)
% em       A positive order to be plotted (random default)
% N        A discretization level (default: 100)
%
% OUTPUT:
%
% errs     A matrix with the errors for everything tested
% l, m     The degree and orders that have been tested
% x        The ordinates at which this is being evaluated
%
% EXAMPLE:
%
% ilk('demo') 
%
% Tests a randomly picked lemma and makes two pictures, one with a
% randomly chosen order for a randomly picked degree, and one with a map
% of the error as a function of the abscissa and the allowable order. 
%
% SEE ALSO: 
%
% PLM, XLM, YLM, LEGENDREDIFF, LIBBRECHT, PAUL
%
% Last modified by plattner-at-princeton.edu, 05/16/2010
% Last modified by fjsimons-at-alum.mit.edu, 02/04/2013

% Which of the lemmata?
defval('lemma',1)

if ~isstr(lemma)

  % The argument at which we are evaluating
  defval('N',500)
  x=linspace(-1,1,N);

  % The degree that we shall be showing
  Lmax=100;
  defval('l',max(round(rand*Lmax),1+[lemma==2]));

  % Let's plot one particular one as far as it's allowed
  defval('em',max(round(rand*(l-1-[lemma==2])),0));

  % Test only the positive orders, which includes potenial recursion to -1
  if em<0; error('Test only orders m>=0'); end

  % All orders
  M=0:l;

  % The Legendre functions and derivatives with respect to acos(x),
  % without (-1)^m, from zero
  [P,dP]=libbrecht(l,x,'sch');

  % Supply the only negative order m=-1, increase array size by one
  P=[(-1)*P(2,:) ;  P];
  dP=[nan(1,N)   ; dP];
  arinc=1;

  % Some of the orders for which we can readily verify the formula
  % that is, not at l==m since the recursion calls for m+1
  m=M(2-arinc:end-1)';

  % Plotting stuff
  epssc=10;
  fs=12;

  switch lemma
   case 1
    % ILK LEMMA 1 quoted by Eshagh (2009a) on page 137
    % \frac{dP_l^m(\cos\theta)}{d\theta}=a_{lm}^1P_l^{m-1}+a_{lm}^2P_l^{m+1},
    % a_{lm}^1=\frac{(l+m)(l-m+1)}{2},\quad a_{lm}^2=-\frac{1}{2}
    % where P_l^m is exactly as Dahlen and Tromp B.71 without (-1)^m
    % We can get this by combining DT B.51 and B.55. Note that GR include
    % the (-1)^m in their definition of Plm, which changes some signs if
    % the orders are one apart. Note that we verify this for Schmidt instead.
    stronk='derivative of Legendre function (Ilk lemma 1)';
    
    % The next two directly from Ilk's First Lemma
    a1=repmat((l+m).*(l-m+1),1,N)/2;
    a2=-1/2;

    % The next three to make the Lemma applicable to Matlab's Schmidt harmonics
    a3=repmat(   sqrt(l+m)./sqrt(l-m  )./sqrt(2-[m  ==0]),1,N);
    a4=repmat(1./sqrt(l-m)./sqrt(l-m+1)./sqrt(2-[m-1==0]),1,N);
    a5=repmat(   sqrt(l+m).*sqrt(l+m+1)./sqrt(2-[m+1==0]),1,N); 

    % This is to be tested for all applicable orders
    left=     a3.*dP(2:end-1,:);
    right=a4.*a1.* P(1:end-2,:)+a5.*a2.*P(3:end,:);
   case 2
    % ILK LEMMA 2 quoted by Eshagh (2009a) on page 137
    % m\frac{P_l^m(\cos\theta)}{d\sin\theta}=
    %                    b_{lm}^1P_{l-1}^{m-1}+b_{lm}^2P_{l-1}^{m+1}, 
    % b_{lm}^1=\frac{(l+m)(l+m-1)}{2},\quad b_{lm}^2=\frac{1}{2}
    % where P_l^m is exactly as Dahlen and Tromp B.71 without (-1)^m
    stronk='m/sin(th) x Legendre function (Ilk lemma 2)';

    % In this Lemma we are comparing to P_{l-1}. Therefore it needs to be
    % calculated too and accordingly normalized
    Plm1=libbrecht(l-1,x,'sch');

    % Some of the orders for which we can readily verify the formula:
    % also not at l==m-1 since the recursion calls for m+1 at l-1
    m=m(1:end-1);
    
    % Supply the only negative order m=-1, increase array size by one
    Plm1=[(-1)*Plm1(2,:) ;  Plm1];
    
    % The next two from the Ilk Lemma
    b1=repmat([m~=0],1,N).*repmat((l+m).*(l+m-1),1,N)/2;
    b2=repmat([m~=0],1,N).*1/2;

    % The next two to make the Lemma applicable to Matlab's Schmidt harmonics
    b3=repmat(m./sqrt(2-[m  ==0]),1,N);
    b4=repmat(1./sqrt(l+m-1)./sqrt(l+m)./sqrt(2-[m-1==0]),1,N);
    b5=repmat(   sqrt(l-m-1).*sqrt(l-m)./sqrt(2-[m+1==0]),1,N);
    
    % P is a function of \cos(x). In the formulae we divide by \sin(x).
    div_sinx=repmat(1./sin(acos(x)),length(m),1);

    % This is to be tested for all applicable orders
    left=b3.*div_sinx.*P   (2:end-2,:);
    right=     b4.*b1.*Plm1(1:end-2,:)+b5.*b2.*Plm1(3:end,:);

    try
      % End points according to L'Hopital's rule 
      left(:,x==-1)=-b3(:,1).*dP(2:end,x==-1);
    end
    try
      left(:,x== 1)= b3(:,1).*dP(2:end-2,x== 1);
    end
   case 3
    % ILK LEMMA 3 quoted by Eshagh (2009a) on page 137
    % m\frac{P_l^m(\cos\theta)}{d\sin\theta}=
    %                 b'_{lm}^1P_{l+1}^{m-1}+b'_{lm}^2P_{l+1}^{m+1},
    % b'_{lm}^1=\frac{(l-m+1)(l-m+2)}{2},\quad b'_{lm}^2=\frac{1}{2}
    % where P_l^m is exactly as Dahlen and Tromp B.71 without (-1)^m
    stronk='m/sin(th) x Legendre function (Ilk lemma 3)';

    % In this Lemma we are comparing to P_{l+1}. Therefore it needs to be
    % calculated too and accordingly normalized
    Plp1=libbrecht(l+1,x,'sch');
    
    % Supply the only negative order m=-1, increase array size by one
    Plp1=[(-1)*Plp1(2,:) ; Plp1];

    % The next two from the Ilk Lemma
    bp1=repmat([m~=0],1,N).*repmat((l-m+1).*(l-m+2),1,N)/2;
    bp2=repmat([m~=0],1,N).*1/2;
    
    % The next three to make the Lemma applicable to Matlab's Schmidt harmonics
    bp3=repmat(m./sqrt(2-[m  ==0]),1,N);
    bp4=repmat(1./sqrt(l-m+2)./sqrt(l-m+1)./sqrt(2-[m-1==0]),1,N);
    bp5=repmat(   sqrt(l+m+2).*sqrt(l+m+1)./sqrt(2-[m+1==0]),1,N);
    
    % P is a function of \cos(x). In the formulae we divide by \sin(x).
    div_sinx=repmat(1./sin(acos(x)),length(m),1);
    
    % This is to be tested for all applicable orders
    left =bp3.*div_sinx.*P   (2:end-1,:);
    right=     bp4.*bp1.*Plp1(1:end-3,:)+bp5.*bp2.*Plp1(3:end-1,:);
    
    try
      % End points according to L'Hopital's rule
      left(:,x==-1)=-bp3(:,1).*dP(2:end-1,x==-1);
    end
    try
      left(:,x== 1)= bp3(:,1).*dP(2:end-1,x== 1);
    end
  end

  % Calculate the error
  errs=left-right;
  err=left(em+arinc,:)-right(em+arinc,:);

  clf
  ah(1)=subplot(211);
  plot(x,left(em+arinc,:),'+')
  hold on
  plot(x,right(em+arinc,:),'r')
  title(sprintf(...
      '%s ; l = %i ; m= %i',stronk,l,em),...
        'FontS',fs)
  hold off

  ah(2)=subplot(212);
  plot(x,err,'+')
  ylim([min(-epssc*eps,-max(abs(err))) max(epssc*eps,max(abs(err)))]);
  title(sprintf(...
      'rms error %s = %6.3g','\epsilon',sqrt(mean(err.^2))),...
        'FontS',fs) 

  nolabels(ah,1)

  % Halt if it is no good; be lenient
  tolex=-round(log10(1e4*epssc*eps));
  difer(err,tolex,[],NaN)
elseif strcmp(lemma,'demo')
  figure(1); clf
  lems=ceil(rand*3);
  [errs,l,m,x]=ilk(lems);
  figure (2); clf; 
  mim=minmax(abs(errs));
  k=kelicol; k(1,:)=[1 1 1];
  imagefnan([x(1) m(1)],[x(end) m(end)],abs(errs),k,mim,[],[],0);
  axis normal
  axis ij; ylabel('order m') ; xlabel('abscissa x'); longticks(gca,2); 
  title(sprintf('lemma % i degree l = %i',lems,l));
  [cb,xcb]=addcb('hor',mim,mim,k,mim(end));
  movev(cb,-0.05)
  set(xcb,'string',sprintf('abs(err), machine eps = %12.8e',eps))
  fig2print(gcf,'portrait')
end

% Optional output
varns={errs,l,m,x};
varargout=varns(1:nargout);

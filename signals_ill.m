function signals_ill(inp)
% Illustrates different signals on a convenient plot

defval('inp',1)

if inp==1
  pin=[1 11 2 3 5 6];
  pon=15;
else
  pin=[8 12 9 10 7 15];
  pon=12;
end  

figure(2)
ah=krijetem(subnum(2,3));

for ind=1:6
  [p(ind),xl(ind),yl(ind),tl(ind)]=plotit(pin,ind,ah);
end

fig2print(gcf,'landscape')
set(tl,'FontS',pon)
longticks(ah)
figdisp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,xl,yl,tl]=plotit(pin,pindex,ah)
figure(1)
[B,F,T,t,wlen,Fs]=signals(pin(pindex),[],'/home/fjsimons/MERMAID/SIGNALS/');
axes(ah(pindex));
p=imagesc(T,F,B);axis xy; colormap(jet)    
xl=xlabel(sprintf('%s ; %3.1f s window','Time (s)',wlen/Fs));
yl=ylabel(sprintf('%s ; Nyquist %5.1f Hz','Frequency (Hz)',Fs/2));
tl=title(t);

function rint=findrad(rad,rs,mmes)
% rint=FINDRAD(rad,rs,mmes)
%
% Nearest-neighbor interpolation of a sorted sequence with repeats.
% If requested point is ON the 'discontinuity', returns the upper side.
%
% INPUT:
%
% rad      A set of radii, sorted, with possible repeats
% rs       A requested radius
% mmes     The type of message to be displayed upon success [default: 1]
%          If 3, then no message at all
%
% OUTPUT:
%
% rint    The index such that rad(rint) is close to rs 
%
% Last modified by fjsimons-at-alum.mit.edu, 08/08/2008

defval('mmes',1)

% Check input was sorted
difer(rad-sort(rad),[],[],NaN)

% Find upper and lower ranges
rl=indeks(find(rad<=rs),'end');
ru=find(rad>rs);
if ~isempty(ru); ru=indeks(ru,1); end

% Take nearest neighbor or upperside and output
[ra,re]=min(abs(rs-rad([ru rl])));
rint=indeks([ru rl],re);

% Diagnostics
switch mmes
 case 3
  return
  % It knows it's a source, and that the input is in meters, and what the
  % value is of an unacceptable distance.
 case 1 
  strmes='Source requested at %7.3f km depth, evaluated at %7.3f depth';
  badval=5000;
 case 2
  strmes='Target requested at %7.3f km depth, evaluated at %7.3f depth';
  badval=5000;
 otherwise
  error('Specify valid message string code')
end

% Protect when it's crazy
disp(sprintf(strmes,(rad(end)-rs)/1000,(rad(end)-rad(rint))/1000))
if ra>badval
  warning('This is more than %7.3g km off the target',ra/1000)
end



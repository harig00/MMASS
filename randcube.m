function X=randcube(d1,d2,d3,fnX)
% X=randcube(d1,d2,d3,fnX)
% 
% Creates a standard cubed-sphere full of rand's
%
% INPUT:
%
% d1,d2,d3   The lengths in the respective dimensions
% fnX        If given, takes the field names from here
%
% OUTPUT:
%
% X          The output structure with six three-dimensional arrays
% 
% SEE ALSO:
%
% CELLNAN, CUBEMATS, NANCUBE
%
% Last modified by fjsimons-at-alum.mit.edu, 02/18/2010

defval('d1',128)
defval('d2',128)
defval('d3',37)
defval('fnX',{'xplus','zminus','yplus','xminus','zplus','yminus'});

for index=1:length(fnX)
  % Uses dynamic structure field variables
  X.(fnX{index})=rand(d1,d2,d3);
end

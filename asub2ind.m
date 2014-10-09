function celnr=asub2ind(lat,dlon,nmr,rowcol)
% celnr=ASUB2IND(lat,dlon,nmr,[rownr colnr])
%
% Matches a equal-area grid row and 'column' number to a sequential cell index
%
% lat         latitudes describing the grid, column vector [degrees]
% dlon        longitudes describing the grid, column vector [degrees]
% nmr         the number of cells per row in the grid
% rownr       the row numbers of the requested cells
% colnr       the 'column' numbers of the requested cells
%
% OUTPUT:
%
% celnr       the requested running cell index vector
%
% SEE ALSO:
%
% AIND2SUB, ACOR2IND
% 
% Last modified by fjsimons-at-alum.mit.edu, 06/13/2007

[rownr,colnr]=deal(rowcol(:,1),rowcol(:,2));
nmr=cumsum([0 ; nmr(:)]);
celnr=nmr(rownr)+colnr;


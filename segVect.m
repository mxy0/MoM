%{
Returns the vector for segment n, given starting and end offsets 

Preconditions:
N: how many segments the antenna is divided into
n: which point index [1, N+1] to grab coordinates 
offset: fraction of segment from point n [-1, 1]
isUnit: 0=return vector, 1=return unit vector

Postconditions:
vect: vector or unit vector of point n
%}
function vect = segVect(N, n, startOffset, endOffset, isUnit)
	[xyzNeg] = segLoc(N, n, startOffset);
	[xyzPos] = segLoc(N, n, endOffset);	
    vect = xyzPos - xyzNeg;
    
    % requesting unit vector
    if isUnit==1
        vect = vect./sqrt(dot(vect,vect));
    end
end

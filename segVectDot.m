%{
returns the electric field dot product with segment unit vector

Preconditions:
N: how many segments the antenna is divided into
n: which point index [1, N+1] to grab coordinates 
startOffset: starting location of vector at point n
endOffset: starting location of vector at point n
incidentVect: vector value of impinging E field

Postconditions:
prod: dot product with segment unit vector and electric field
%}
function prod = segVectDot(N, m, startOffset, endOffset, incidentVect)
    sampleSize = length(startOffset);
    prod = zeros(sampleSize,1);
    
    for i = 1:sampleSize
%        vect = segVect(N, m, startOffset(i), endOffset(i), 1); 
       vect = segVect(N, m, startOffset(i), endOffset(i), 1);    
       prod(i) = dot(vect, incidentVect);
    end
end
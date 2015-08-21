%{
Returns the distance between the points m and n 

Preconditions:
N: how many segments the antenna is divided into
m: which point index [1, N+1] to grab coordinates of source charge
offsetM: offset from point m along the segment [-1,1]
n: which point index [1, N+1] to grab coordinates of observation point
offsetN: offset from point n along the segment [-1,1]

Postconditions:
RDist: scalar distance between points m and n including wire radius
%}
function [RDist] = Rmn(N, m, offsetM, n, offsetN)                      
	[xyzM] = segLoc(N, m, offsetM);
	[xyzN, a] = segLoc(N, n, offsetN);
	RDist = sqrt((a)^2+sum((xyzM-xyzN).^2));
end
%{
Returns coordinates of segment n
Preconditions:
N: how many segments the antenna is divided into
n: which point index [1, N+1] to grab coordinates 
offset: fraction of segment from point n [-1, 1]

Postconditions:
xyz: [x y z] coordinate of the point
a: wire radius of structure

NOTE: change for different wire antenna geometries
%}
function [xyz, a] = segLoc(N, n, offset)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add following line in all code!!!!
    nSel = n-1+offset;      % which point of segment to return coordinates
                            % n-1 gets to point location

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % antenna design parameters
    h = 1e-6;         	% total antenna length
    gamma1 = h/2;
    gamma2 = h/2;
    L = gamma1+gamma2;
    psi_ang = 90*2*pi/360;  % opening angle of antenna
    a = 0.001e-6;             % thin wire approximation radial width
                            % 200 nm wide and 100 nm thick antenna

%     assignin('base', 'L', L);   % expose variable value to workspace    
%     assignin('base', 'a', a);   % expose variable value to workspace
    
    corner = N/2;  
    
    x = 0;
    y = L/N*(nSel-corner);
    z = 0;
    
    xyz = [x y z];
end

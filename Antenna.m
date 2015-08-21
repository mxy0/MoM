classdef Antenna
    properties (Access = public)
        mu_0 = 4*pi*1e-7;       % vacuum permeability
        epi_0 = 8.854*1e-12;    % permittivity
        c                       % speed of light    
        
        % antenna division
        N                       % how many segments to divide antenna into
        M = 6                   % number of sample points to take for triangle basis 
                                % must be even
                                
        mu_r                    % relative permeability
        epi_r                   % relative permittivity
        k                       % wavenumber
        omega                   % operating frequency
        E_inc                   % incident field magnitude
        E_vect                  % incident field vector
        
        Z                       % full impedance matrix for MoM
        V                       % potential matrix along wire
        I                       % current distribution along wire
        ZWire                   % impedance per segment along wire
        inputIm                 % input impedance
        ind                     % antenna input inductance
        QFactor                 % q factor 
        del_l                   % segment length 
        a                       % wire radius
        
        assump                  % assume end currents are zero
    end
    
    methods
        % constructor
        function obj = Antenna(dir, N, freq, epi_r, mu_r, E_inc, E_vect, assump)
            added_path = [pwd '\' dir];           % choose folder for antenna geometry
            addpath(added_path);
            
            obj.N = N;
            obj.c =  1/sqrt(obj.mu_0*obj.epi_0);  % vacuum speed of light
            obj.epi_r = epi_r;
            obj.mu_r = mu_r;
            obj.omega = 2*pi*freq;                % operating frequency
            obj.k = obj.omega*sqrt(obj.epi_0*epi_r*obj.mu_0*mu_r); % k vector
            
            obj.E_inc = E_inc;
            obj.E_vect = E_vect;
            
            obj.assump = assump;
            
            % length of each segment
            del_l = segVect(obj.N, 1, 0, 1, 0);        
            obj.del_l = sqrt(dot(del_l, del_l));                % magnitude
            % retreive wire radius
            [~, obj.a] = segLoc(obj.N, 1, 0);
            obj.Z = calcZ(obj);
            obj.V = calcV(obj);
            obj.I = calcI(obj, assump);
           
            obj.ZWire = calcZWire(obj);         % impedance of antenna, per segment
            obj.inputIm = calcInputIm(obj, assump);
            obj.QFactor = calcQFactor(obj);     % Q factor of antenna
            obj.ind = calcInd(obj);             % inductance of antenna (j*omega*L)

        end
        
        % calculate the MoM impedance matrix
        function Z = calcZ(obj)
            % copy object's properties for readability
            N = obj.N;
            M = obj.M;
            del_l = obj.del_l;
            k = obj.k;
            a = obj.a;
            omega = obj.omega;
            epi_0 = obj.epi_0;
            mu_0 = obj.mu_0;
            
            % matrices to track row and column indices of each element for basis and
            % testing functions
            [colCountM, rowCountM] = meshgrid(1:M, 1:M);
            % normalize column and row matrixes so that [1,M] -> [-1+1/M,1-1/M]
            % needed to use for integral approximation using riemann sums
            colCountMNorm = norm(colCountM, M);
            rowCountMNorm = norm(rowCountM, M);

            basisProd = arrayfun(@f, rowCountMNorm).*arrayfun(@f, colCountMNorm);

            % derivative of triangle basis and testing functions
            % p is [-1,1]
            fDer = @(p) sign(-p)*1/del_l;
            basisDerProd = arrayfun(fDer, rowCountMNorm).*arrayfun(fDer, colCountMNorm);
            
            % m = observation point
            % n = source point
            % calculate Zmn impedance
            
            % m != n Kernal
            K_pq = @(m, offsetM, n, offsetN) exp(-1i*k*Rmn(N, m, offsetM, n, offsetN))...
                    /(4*pi*Rmn(N, m, offsetM, n, offsetN));
            G_pq = @(m,n) arrayfun(K_pq, repmat(m,M), rowCountMNorm, repmat(n,M), colCountMNorm);


            % segment unit vector dot product
            segUnitDot = @(m, offsetM, n, offsetN) dot(segVect(N, m, offsetM-1/(M+1), offsetM+1/(M+1), 1), ...
                segVect(N, n, offsetN-1/(M+1), offsetN+1/(M+1), 1));
            dotProd = @(m,n) arrayfun(segUnitDot, repmat(m,M), rowCountMNorm,...
                repmat(n,M), colCountMNorm);

            % width of each Riemann sum column
            sampleWidth = repmat(2*del_l/M, M);
            % first and last segments require half triangles as bases
            startEdge = [zeros(M/2,M); ones(M/2, M)];
            endEdge = [ones(M/2, M); zeros(M/2, M)];
            
            % impedance MxM matrix terms
            ZmnMatrix = @(m,n) G_pq(m,n).*(sampleWidth.^2).*(1i*omega*mu_0*basisProd.*dotProd(m,n)...
                -1i/(omega*epi_0)*basisDerProd);
            
            % impedance NxN matrix terms
            Zmn = @(m,n) sum(sum(ZmnMatrix(m,n)));
            % special scenerios for beginning and end segments due to 
            % half triangle approximation
            ZmnTop = @(m,n) sum(sum(startEdge.*ZmnMatrix(m,n)));
            ZmnBottom = @(m,n) sum(sum(endEdge.*ZmnMatrix(m,n)));
            ZmnLeft = @(m,n) sum(sum(transpose(startEdge).*ZmnMatrix(m,n)));
            ZmnRight = @(m,n) sum(sum(transpose(endEdge).*ZmnMatrix(m,n)));
            ZmnTopLeft = @(m,n) sum(sum(startEdge.*transpose(startEdge).*ZmnMatrix(m,n)));
            ZmnTopRight = @(m,n) sum(sum(startEdge.*transpose(endEdge).*ZmnMatrix(m,n)));
            ZmnBottomLeft = @(m,n) sum(sum(endEdge.*transpose(startEdge).*ZmnMatrix(m,n)));
            ZmnBottomRight = @(m,n) sum(sum(endEdge.*transpose(endEdge).*ZmnMatrix(m,n)));

            % matrices to track row and column indices of each basis function
            [colCountN, rowCountN] = meshgrid(2:N, 2:N);

            Z = zeros(N+1,N+1);
            % unmodified Z matrix with no self terms
            Z(2:N, 2:N) = arrayfun(Zmn, rowCountN, colCountN);
            Z(1, 2:N) = arrayfun(ZmnTop, repmat(1, 1, N-1), 2:N);
            Z(N+1, 2:N) = arrayfun(ZmnBottom, repmat(N+1, 1, N-1), 2:N);
            Z(2:N, 1) = arrayfun(ZmnLeft, 2:N, repmat(1, 1, N-1));
            Z(2:N, N+1) = arrayfun(ZmnRight, 2:N, repmat(N+1, 1, N-1));
            Z(1,1) = ZmnTopLeft(1,1);
            Z(1, N+1) = ZmnTopRight(1, N+1);
            Z(N+1, 1) = ZmnBottomLeft(N+1, 1);
            Z(N+1, N+1) = ZmnBottomRight(N+1, N+1);
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % calculate self term contributions
            % basis for an ascending triangle [0,1]
            basisAsc =  norm([M/2+1:M],M);
            sampleWidth1 = repmat(2*del_l/M, 1, length(basisAsc));
            % ascending triangle function [0,del_l]
            selfCount = basisAsc*del_l;

            % equivalent radius for self-terms
            a=a*(1-0.40976*a/del_l);
            
            % x is observation point[-del_l, del_l]
            S1 = @(x) 1/del_l*sqrt(a^2+(x-del_l)^2)-1/del_l*sqrt(a^2+x^2)+...
                x/del_l*log((sqrt(a^2+(x-del_l)^2)-x+del_l)/(sqrt(a^2+x^2)-x))-1i*k*del_l/2;

%             S1 = @(x) 1/del_l*sqrt(a^2+(x-del_l)^2)-1/del_l*sqrt(a^2+x^2)+...
%                 x/del_l*log((sqrt(a^2+x^2)+x)/(sqrt(a^2+(x-del_l)^2)+x-del_l))-1i*k*del_l/2;
            
            % due to symmetry, we care only about the ascending triangle 
            selfS1 = arrayfun(S1, selfCount);

            S2 = @(x) 1/(del_l^2)*(log((x+sqrt(a^2+x^2))/(x-del_l+sqrt(a^2+(x-del_l)^2)))-1i*k*del_l);
            selfS2 = arrayfun(S2, selfCount);

            % -1 slope when derivatives of source and testing functions are of
            % opposing sign
            % 1 slope when derivatives of source and testing functions are of
            % the same sign
            selfMMatrix1 = 1/(4*pi)*sampleWidth1.*...
                (1i*omega*mu_0*arrayfun(@f, basisAsc).*selfS1-1i/(omega*epi_0)*(-selfS2));

            selfMMatrix2 = 1/(4*pi)*sampleWidth1.*...
                (1i*omega*mu_0*arrayfun(@f, -(1-basisAsc)).*selfS1-1i/(omega*epi_0)*(selfS2));

            % m = observation point
            % n = source point
            % calculate Zmn impedance
            
            % which terms in the NxN matrix have overlapping segments
            mnOverlapN = full(spdiags(ones(N+1), 1, N+1, N+1));       % m+1=n
            nnOverlapN = full(spdiags(ones(N+1), 0, N+1, N+1));       % m=n
            nmOverlapN = full(spdiags(ones(N+1), -1, N+1, N+1));       % m=n+1

            % which terms in the MxM sampling matrix correspond to overlapping segments
            % all observation points for given source segment
            mnOverlapM = [zeros(M/2, M); ones(M/2, M/2), zeros(M/2, M/2)];
            nnOverlapM = [ones(M/2, M/2), zeros(M/2,M/2); zeros(M/2,M/2), ones(M/2, M/2)];
            nmOverlapM = [zeros(M/2, M/2), ones(M/2, M/2); zeros(M/2, M)];

            % calculate Z terms with overlapping segments
            ZmnOverlap = @(m,n) sum(sum(ZmnMatrix(m,n).*(1-mnOverlapM)))+sum(sum(selfMMatrix1));
            ZnnOverlap = @(m,n) sum(sum(ZmnMatrix(m,n).*(1-nnOverlapM)))+2*sum(sum(selfMMatrix2));
            ZnmOverlap = @(m,n) sum(sum(ZmnMatrix(m,n).*(1-nmOverlapM)))+sum(sum(selfMMatrix1));

            ZOverlap = zeros(N+1);
            % recalculate all matrix values that have self terms
            for i = 2:N
                ZOverlap(i, i-1) = ZnmOverlap(i, i-1);
                ZOverlap(i, i+1) = ZmnOverlap(i, i+1);
                ZOverlap(i, i) = ZnnOverlap(i, i);   
            end
            % beginning and end segments with the half triangle basis must
            % be specially calculated
            ZOverlap(1,1) = sum(sum(selfMMatrix2));
            ZOverlap(1,2) = sum(sum(ZmnMatrix(1,2).*(1-mnOverlapM).*startEdge))+sum(sum(selfMMatrix1));
            ZOverlap(2,1) = sum(sum(ZmnMatrix(2,1).*(1-nmOverlapM).*transpose(startEdge)))+sum(sum(selfMMatrix1));
            ZOverlap(N+1,N+1) = sum(sum(selfMMatrix2));
            ZOverlap(N,N+1) = sum(sum(ZmnMatrix(N,N+1).*(1-mnOverlapM).*transpose(endEdge)))+sum(sum(selfMMatrix1));
            ZOverlap(N+1,N) = sum(sum(ZmnMatrix(N+1,N).*(1-nmOverlapM).*endEdge))+sum(sum(selfMMatrix1));
            
            assignin('base', 'ZBefore', Z);   % expose variable value to workspace 

            % complete Z matrix with updated self terms
            Z = Z.*(1-(mnOverlapN+nnOverlapN+nmOverlapN))+ ZOverlap;

            ZAbs = abs(Z);

            assignin('base', 'ZAbs', ZAbs);   % expose variable value to workspace 
            assignin('base', 'ZOverlap', ZOverlap);   % expose variable value to workspace 
        end
        
        % calculate the antenna's potential 
        function V = calcV(obj)
            omega = obj.omega;
            mu_0 = obj.mu_0;
            E_inc = obj.E_inc;
            E_vect = obj.E_vect;
            M = obj.M;
            N = obj.N;
            del_l = obj.del_l;
            
            rowCountMNorm = norm(transpose(1:M), M);
            sampleWidth = repmat(2*del_l/M, M,1);
            % potential boundary matrix
            basis = arrayfun(@f, rowCountMNorm);

            bm = @(m) sum(sum(sampleWidth.*(basis.*segVectDot(N, m,...
                rowCountMNorm(:,1)-1/(M+1), rowCountMNorm(:,1)+1/(M+1), E_inc*E_vect))));
            V = arrayfun(bm,transpose(1:N+1));
            V(1) = V(1)/2;
            V(end) = V(end)/2;
            
         
            V = zeros(N+1,1);
            bm1 = sum(sum(sampleWidth.*(basis*1/del_l)));
            
            V(ceil((N+2)/2)) = V(ceil((N+2)/2)) + bm1/2;
            V(floor((N+2)/2)) = V(floor((N+2)/2)) + bm1/2;

%             V(1) = bm1/2;
%             V(end) = bm1/2;

%             V(1) = bm1/2;
%             V(end) = -bm1/2;
            
%             V(1) = bm1;
        end
        
        % calculate current distribution
        function I = calcI(obj, typeAssume)
            N = obj.N;
            
            I = zeros(N+1,1);
            % assume end currents are zero
            if typeAssume == 0
                I(2:N) =  obj.Z(2:N, 2:N)\obj.V(2:N); 
            else
                I =  obj.Z\obj.V; 
            end

        end
        
        % calculate impedance of each segment along antenna
        function ZWire = calcZWire(obj)
           ZWire = obj.V./obj.I; 
        end
        
        % calculate total q factor
        function inputIm = calcInputIm(obj, assump)
            N = obj.N;
            
            if assump == 0
                inputIm = sum(obj.ZWire(2:N)); 
            else
                inputIm = sum(obj.ZWire);
            end
        end
        
        % calculate total q factor
        function QFactor = calcQFactor(obj)
            QFactor = abs(imag(obj.inputIm)/real(obj.inputIm));
        end
        
        % calculate inductance of antenna
        function ind = calcInd(obj)
           ind = -imag(obj.inputIm)/obj.omega;
        end
        
        % plot current distribution
        function plotCurr(obj)
            I = obj.I;
            N = obj.N;
%             I = [0; obj.I; 0];
            %plot current distribution
            figure(1);
%             plot(1/(length(I)-1)*(0:length(I)-1), abs(I));
            plot(((1:length(I))-1/2)/(N+1), abs(I));
            xlabel('Segment Location');
            ylabel('Current (I)');
            % title('\bf blah')
        end

        % plot wire location and geometry
        function plotWire(obj)
            % plot the wire geometry
            N = obj.N;
            E_vect = obj.E_vect;
            x_coor = zeros(N+1,1);
            y_coor = zeros(N+1,1);
            z_coor = zeros(N+1,1);
            
            for i=1:N+1
                coor = segLoc(N, i, 0);
                x_coor(i,1) = coor(1);
                y_coor(i,1) = coor(2);
                z_coor(i,1) = coor(3);
            end

            factor = max(max(max(abs(x_coor)), max(abs(y_coor))), max(abs(z_coor)));
            E_x = [0 E_vect(1)].*factor;
            E_y = [0 E_vect(2)].*factor;
            E_z = [0 E_vect(3)].*factor;

            figure(2);
            plot3(x_coor, y_coor, z_coor);
%             plot3(x_coor, y_coor, z_coor, E_x, E_y, E_z);
            hold on;
%             plot3(E_x(2),E_y(2), E_z(2),'r.','MarkerSize',15)
            xlabel('x');
            ylabel('y');
            zlabel('z');
            axis([-8e-6,8e-6,-8e-6,8e-6,-8e-6,8e-6]);
%             axis([-1.5e-6,1.5e-6,-1.5e-6,1.5e-6,-1.5e-6,1.5e-6]);
            %plot(x_coor, y_coor, E_x, E_y);
            daspect([1 1 1]);
            grid on;
        end

        % plot field impinging on plane
        % fType: 1-calculate h field, 0-calculate e field
        function plotField(obj, xIn, yIn, zIn, fType)
            N = obj.N;
            del_l = obj.del_l;
            omega = obj.omega;
            mu = obj.mu_0*obj.mu_r;
            k = obj.k;
            I = [0; obj.I; 0];
            a = obj.a;
            epi = obj.epi_r*obj.epi_0;
            
            % calculate radiated field components
            distSO = @(m, obsLoc) sqrt(sum((obsLoc - segLoc(N, m, 0)).^2));
            
            % plot E-field
            if fType == 0
                fieldName = 'E';
                % far field 
                FmFar = @(m, obsLoc) -1i*omega*mu/(4*pi)*I(m)/(2*pi*a)*segVect(N, m, -del_l/2, del_l/2, 0)*...
                    exp(-1i*k*distSO(m, obsLoc))/distSO(m, obsLoc);

                % near field
                G1 = @(r) (-1-1i*k*r+k^2*r^2)/(4*pi*r^3);
                G2 = @(r) (3+3*1i*k*r-k^2*r^2)/(4*pi*r^5);
                B = @(m, obsLoc) sum((obsLoc - segLoc(N, m, 0)).*...
                    (I(m)/(2*pi*a)*segVect(N, m, -del_l/2, del_l/2, 0)));
                FmNear = @(m, obsLoc) 1/(1i*omega*epi)*(G1(distSO(m, obsLoc))*...
                    I(m)/(2*pi*a)*segVect(N, m, -del_l/2, del_l/2, 0) + ...
                    B(m, obsLoc)*(obsLoc - segLoc(N, m, 0))*G2(distSO(m, obsLoc)))*...
                    exp(-1i*k*distSO(m, obsLoc));
            % plot H-field
            elseif fType == 1
                fieldName = 'H';
                % far field
                FmFar = @(m, obsLoc) sqrt(mu*epi)/(4*pi)*cross(obsLoc/dot(obsLoc, obsLoc),...
                    I(m)/(2*pi*a)*segVect(N, m, -del_l/2, del_l/2, 0))*...
                    exp(-1i*k*distSO(m, obsLoc))/distSO(m, obsLoc);
                % near field
                FmNear = @(m, obsLoc) -cross(obsLoc - segLoc(N, m, 0), ...
                    I(m)/(2*pi*a)*segVect(N, m, -del_l/2, del_l/2, 0))*...
                    (1+1i*k*distSO(m, obsLoc))/(4*pi*distSO(m, obsLoc)^3)*...
                    exp(-1i*k*distSO(m, obsLoc));
            end
            
            FrFullMatrix = zeros(N,3);

            if length(xIn) == 1
                inLoop = length(yIn);
                outLoop = length(zIn);
                xGrid = repmat(xIn, outLoop, inLoop);
                yGrid = repmat(yIn, outLoop, 1);
                zGrid = repmat(transpose(zIn), 1, inLoop);
            elseif length(yIn) ==1
                inLoop = length(xIn);
                outLoop = length(zIn);
                xGrid = repmat(xIn, outLoop, 1);
                yGrid = repmat(yIn, outLoop, inLoop);
                zGrid = repmat(transpose(zIn), 1, inLoop);
            elseif length(zIn) == 1
                inLoop = length(xIn);
                outLoop = length(yIn);
                xGrid = repmat(xIn, outLoop, 1);
                yGrid = repmat(transpose(yIn), 1, inLoop);
                zGrid = repmat(zIn, outLoop, inLoop);
            end
            assignin('base', 'xGrid', xGrid);
            assignin('base', 'yGrid', yGrid);
            assignin('base', 'zGrid', zGrid);
            FrMag = zeros(outLoop, inLoop);
            FrX = zeros(outLoop, inLoop);
            FrY = zeros(outLoop, inLoop);
            FrZ = zeros(outLoop, inLoop);
            
            for l = 1:outLoop
                for j = 1:inLoop
                	obsLoc = [xGrid(l, j),yGrid(l, j),zGrid(l, j)];
                    for i = 1:N
                        dist = distSO(i+1, obsLoc);
                        if dist <= a
                            FrFullMatrix(:,:) = 0;
                            break;
                        end
                        
                        if k*distSO(i+1, obsLoc) >= 10
                            FrFullMatrix(i,:) = FmFar(i+1, obsLoc);
                        else
                            FrFullMatrix(i,:) = FmNear(i+1, obsLoc);
                        end
                    end
%                     Er = exp(1i*60*pi/180)*sum(ErFullMatrix);
                    Fr = sum(FrFullMatrix);
                    FrReal = real(Fr);
                    FrX(l, j) = FrReal(1);
                    FrY(l, j) = FrReal(2);
                    FrZ(l, j) = FrReal(3);
                end
            end
%             FrMag = sqrt(FrX.^2+FrY.^2+FrZ.^2)*1e3;
            FrMag = FrZ;
            
            assignin('base', 'FrX', FrX);   % expose variable value to workspace 
            assignin('base', 'FrY', FrY);   % expose variable value to workspace 
            assignin('base', 'FrZ', FrZ);   % expose variable value to workspace 
            
            assignin('base', 'FrMag', FrMag);   % expose variable value to workspace 

            figure(3);
            if length(xIn) == 1
                surf(yIn, zIn, FrMag, 'EdgeColor', 'None');
                set(gca,'XLim',[min(yIn) max(yIn)],'YLim',[min(zIn) max(zIn)])
                xlabel('y');
                ylabel('z');
            elseif length(yIn) ==1
                surf(xIn, zIn, FrMag, 'EdgeColor', 'None');
                set(gca,'XLim',[min(xIn) max(xIn)],'YLim',[min(zIn) max(zIn)])
                xlabel('x');
                ylabel('z');
            elseif length(zIn) == 1
                surf(xIn, yIn, FrMag, 'EdgeColor', 'None');
                set(gca,'XLim',[min(xIn) max(xIn)],'YLim',[min(yIn) max(yIn)])
                xlabel('x');
                ylabel('y');
            end

            zlabel([fieldName '-Field Magnitude']);
            aspectRatio = daspect;
            aspectRatio(1) = min(aspectRatio);
            aspectRatio(2) = min(aspectRatio);
            daspect(aspectRatio);
            colormap(b2r(-max(max(abs(FrMag))),max(max(abs(FrMag)))));
            colorbar('eastoutside');
            set(gca,'TickDir','out')
            box on;
        end
    end
end

% calculate amplitude of a triangle basis - p is [-1,1]
function   triangleBasis = f(p)         
    triangleBasis = 1-abs(p);
end

% normalize a matrix for M sampling points
function normMat = norm(matrix, M)
    normMat = matrix*2/M-(M+1)/M;
end



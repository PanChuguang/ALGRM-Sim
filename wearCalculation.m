function newDisBR = wearCalculation(oldDisBR,qBnext,qBcur,eBJnext,eBJcur,qJcur,qJnext,Fn,Lc,Rj,Kw)
    % function for calculation of wearing effect between [t_{n},t_{n+1}]
    %
    % Input arguments:
    %   oldDisBR -- discrete radius of bearing contour before wearing
    %   qBnext -- generalzied coordinates of bearing body at time step n+1
    %   qBcur -- generalized coordinates of bearing body at time step n
    %   eRJnext -- eccentricity vector from bearing to journal (GCS) n+1
    %   eRJcur -- eccentricity vector from bearing to journal (GCS) n
    %   qJcur -- generalized coordinates of journal body at time step n
    %   qJnext -- generalized coordinates of journal body at time step n+1
    %   Fn -- normal contact force in current time step n
    %   Lc -- journal contact length
    %   Rj -- constant radius of journal
    %   Kw -- Archard's coefficient

    % Output arguments:
    %   newDisBR -- updated bearing radius

    arguments (Input)
        oldDisBR (:,1) double {mustBeFinite,mustBeNonNan}
        qBnext (4,1) double {mustBeFinite,mustBeNonNan}
        qBcur (4,1) double {mustBeFinite,mustBeNonNan}
        eBJnext (2,1) double {mustBeFinite,mustBeNonNan}
        eBJcur (2,1) double {mustBeFinite,mustBeNonNan}
        qJcur (4,1) double {mustBeFinite,mustBeNonNan}
        qJnext (4,1) double {mustBeFinite,mustBeNonNan}
        Fn (1,1) double {mustBePositive,mustBeNonNan,mustBeNonNan}
        Lc (1,1) double {mustBePositive,mustBeNonNan,mustBeFinite}
        Rj (1,1) double {mustBePositive,mustBeNonNan,mustBeFinite}
        Kw (1,1) double {mustBePositive,mustBeNonNan,mustBeFinite}
    end

    arguments (Output)
        newDisBR (:,1) double {mustBeFinite,mustBeNonNan}
    end

    % contact pressure calculation of contacted points
    % obtain contact points firstly
    [~,beta_s,beta_a,beta_b,~,~,~] = getContactRegion(oldDisBR,qBcur,qJcur,eBJcur,Rj);

    % get max contact pressure p_max
    %
    N = numel(oldDisBR); % number of discrete points along the bearing contour
    betas = 0:2*pi/N:2*pi-2*pi/N;
    Rb = @(x) spline([betas,2*pi],[0;oldDisBR;oldDisBR(1);0],mod(x,2*pi));
    segARange = unwrap([beta_a beta_s]);
    segBRange = unwrap([beta_s beta_b]);
    % segASamples = linspace(segARange(1),segARange(2),50);
    % segBSamples = linspace(segBRange(1),segBRange(2),50);
    % arcA = trapz(segASamples,Rb(segASamples));
    % arcB = trapz(segBSamples,Rb(segBSamples));
    arcA = integral(Rb,segARange(1),segARange(2));
    arcB = integral(Rb,segBRange(1),segBRange(2));

    pmax = 4.*Fn./((arcA + arcB)*pi*Lc);
    
    % obatin slip distance from the current time step n
    [~,beta_s_new,~,~,~,~,~] = getContactRegion(oldDisBR,qBnext,qJnext,eBJnext,Rj);
    deltaS = abs(diff(unwrap([beta_s beta_s_new])))*Rb(beta_s_new); % slip distance
   
    % contact points wear calculation
    newDisBR = oldDisBR; % preallocation output variable

    for idx = 1:N
        if isAngleInRange(betas(idx),beta_a,beta_s)
            % contact point is in the front segment
            pk = pmax*sqrt(1-diff(unwrap([betas(idx) beta_s])).^2./diff(unwrap([beta_a beta_s])).^2);
            newDisBR(idx) = newDisBR(idx) + Kw*pk*deltaS;
        elseif isAngleInRange(betas(idx),beta_s,beta_b)
            % contact point is in the back segment
            pk = pmax*sqrt(1-diff(unwrap([betas(idx) beta_s])).^2./diff(unwrap([beta_s beta_b])).^2);
            newDisBR(idx) = newDisBR(idx) + Kw*pk*deltaS;
        end
    end

end




function result = isAngleInRange(x,lb,ub)
    % Determine whether the angle x is located in the [lb, ub] interval
    % The interval may exist angle wrap, e.g., deg2rad([350,10])
    if lb <= ub
        result = (x >= lb) && (x <= ub);
    else
        result = (x >= lb) || (x <= ub);
    end
end
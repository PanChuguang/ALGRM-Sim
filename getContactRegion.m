function  [delta_c,beta_c,beta_s,beta_f,loc_B,loc_P,isContact] = getContactRegion(disBR,qB,qJ,eBJ,Rj)
    % function for obtaining contact region on nonuniform bearing contour
    % 
    % Input arguments:
    %   disBR -> discrete radius of bearing contour
    %   qB -> generalized coordinates of bearing body
    %   qJ -> generalized coordinates of journal body
    %   eBJ -> eccentricity vector from bearing to journal (GCS)
    %   Rj -> constant radius of journal
    %
    % Output arguments:
    %   delta_c -> penetration depth in contact angle \beta_c
    %   beta_c -> contact angle of contact point
    %   beta_s -> entry angle of contact region
    %   beta_f -> exit angle of contact region Note: beta_s may be > beta_f
    %   loc_B -> local coordinates of the contact point on the bearing
    %   loc_P -> local coordinates of the contact point on the journal
    %   isContatc -> contact flag

    arguments (Input)
        disBR (:,1) double {mustBeFinite,mustBeNonNan}
        qB (4,1) double {mustBeFinite,mustBeNonNan}
        qJ (4,1) double {mustBeFinite,mustBeNonNan}
        eBJ (2,1) double {mustBeFinite,mustBeNonNan} 
        Rj (1,1) double {mustBeFinite,mustBeNonNan,mustBePositive}
    end

    arguments (Output)
        delta_c (1,1) double {mustBeFinite,mustBeNonNan}
        beta_c  (1,1) double {mustBeFinite,mustBeNonNan}
        beta_s  (1,1) double {mustBePositive,mustBeFinite}
        beta_f  (1,1) double {mustBePositive,mustBeFinite}
        loc_B  (2,1) double {mustBeFinite,mustBeNonNan}
        loc_P  (2,1) double {mustBeFinite,mustBeNonNan}
        isContact  logical {mustBeScalarOrEmpty,mustBeNonempty}
    end
    
    N = numel(disBR); % number of discretization points on bearing contour
    dbeta = 2*pi/N; % delta angle between two adjacent discrete point
    betas = (0:dbeta:2*pi-dbeta).'; % uniform distributed angles of discrete point
    Rb = @(ang) spline([betas;2*pi],[0;disBR;disBR(1);0],mod(ang,2*pi));

  
    sB = @(ang) consTranMat2D(Rb(ang).*[cos(ang);sin(ang)],"vector")*qB;
    dB = @(ang) sB(ang) - eBJ;
    penDepth = @(ang) Rj - vecnorm(dB(ang));
    
    % Contact condition determine
    isContact = false;
    for idx = 1:N
        if penDepth(betas(idx)) > 0
            isContact = true;
            break;
        end
    end
    
    if ~isContact
        delta_c = nan;
        beta_c = nan;
        beta_s = nan;
        beta_f = nan;
        loc_B = nan;
        loc_P = nan;
        return;
    end

    % obtain contact parameters
    disPD = zeros(N+1,1);

    augBetas = linspace(0,2*pi,N+1);
     for idx = 1:N+1
         disPD(idx) = penDepth(augBetas(idx));
     end

    beta_s0idx = find(diff(sign(disPD))>0,1,"first");
    beta_f0idx = find(diff(sign(disPD))<0,1,"first");
    beta_s1idx = beta_s0idx + 1;
    beta_f1idx = beta_f0idx + 1;
    
    beta_s = fzero(penDepth,[augBetas(beta_s0idx),augBetas(beta_s1idx)]);
    beta_f = fzero(penDepth,[augBetas(beta_f0idx),augBetas(beta_f1idx)]);

    % calculation of contact angle \beta_c
    unwrapRange = unwrap([beta_s;beta_f]);
    intSamples = linspace(unwrapRange(1),unwrapRange(2)); 

    intTargets = zeros(numel(intSamples),1);
    for sampleIdx = 1:numel(intTargets)
        intTargets(sampleIdx) = penDepth(intSamples(sampleIdx));
    end
    
    beta_c = trapz(intSamples,intSamples.*intTargets)./trapz(intSamples,intTargets);
    beta_c = mod(beta_c,2*pi); % unwrap angle to [0,2*pi]
    delta_c = penDepth(beta_c);

    loc_B = Rb(beta_c)*[cos(beta_c);sin(beta_c)];
    loc_P = [dot(dB(beta_c),normalize(qJ(3:4),"norm"));...
            dot(dB(beta_c),[0 -1;1 0]*normalize(qJ(3:4),"norm"))];

end
function [F,locB,locP,Fn] = getContactForce(disBR,qB,dqB,eBJ,deBJ,qJ,predeldot,Rj,K,ce,muf,nu0,nu1)
    % calculation of contact force between bearing and journal with 
    % action point at contact point
    %
    % Input arguments:
    %   disBR -- discrete radius of bearing contour
    %   qB -- generalized coordinates of bearing body
    %   dqB -- generalized coordinates' dot of bearing body
    %   eBJ -- eccentricity vector from bearing to journal (GCS)
    %   deBJ -- derivative of eBJ with respect to time
    %   predeldot -- initial impact velocity (from free flight to contact)
    %   Rj -- constant radius of journal
    %   K -- stiffness of the Hertz model
    %   ce -- coeffcient of restitution 
    %   muf -- friction coefficient
    %   nu0 -- friction tolerance lower bound
    %   nu1 -- friction tolerance upper bound
    % Output arguments:
    %   F -- contact force vector from journal to bearing
    %   locB -- localized contact point coordinates of bearing 
    %   locP -- localized contact point coordinates of journal
    %   Fn -- magnitude of normal contact force

    arguments (Input)
        disBR (:,1) double {mustBeFinite,mustBeNonNan}
        qB (4,1) double {mustBeFinite,mustBeNonNan}
        dqB (4,1) double {mustBeFinite,mustBeNonNan}
        eBJ (2,1) double {mustBeFinite,mustBeNonNan}
        deBJ (2,1) double {mustBeFinite,mustBeNonNan}
        qJ (4,1) double {mustBeFinite,mustBeNonNan}
        predeldot (1,1) double {mustBeFinite,mustBeNonNan}
        Rj (1,1) double {mustBePositive,mustBeFinite,mustBeNonNan}
        K (1,1) double {mustBePositive,mustBeNonNan}
        ce (1,1) double {mustBePositive,mustBeNonNan}
        muf (1,1) double {mustBeFinite,mustBeNonNan,mustBePositive}
        nu0 (1,1) double {mustBeFinite,mustBeNonNan,mustBePositive}
        nu1 (1,1) double {mustBeFinite,mustBeNonNan,mustBePositive}
    end

    arguments (Output)
        F (2,1) double {mustBeFinite,mustBeNonNan}
        locB (2,1) double {mustBeFinite,mustBeNonNan}
        locP (2,1) double {mustBeFinite,mustBeNonNan}
        Fn  (1,1) double {mustBeNonnegative,mustBeFinite}
    end

    % get contact point
    [deltaS,betaS,~,~,locB,locP,~] = getContactRegion(disBR,qB,qJ,eBJ,Rj);

    if isnan(deltaS)
        [F,locB,locP] = deal(zeros(2,1));
        Fn = 0;
        return;
    end

    N = numel(disBR);
    angles = 0:2*pi/N:2*pi-2*pi/N;
    Rb = @(x) spline([angles,2*pi],[0;disBR;disBR(1);0],mod(x,2*pi));
    dvec = consTranMat2D(Rb(betaS).*[cos(betaS);sin(betaS)],"vector")*qB - eBJ;
    dvecdot = consTranMat2D(Rb(betaS).*[cos(betaS);sin(betaS)],"vector")*dqB - deBJ;

    ddeltaS = -dot(normalize(dvec,"norm"),dvecdot);

    % normal force calculation
    Fn = K.*max(0,deltaS).^1.5.*(1+3*(1-ce^2)*ddeltaS/(4*predeldot));
    normvec = normalize(dvec,"norm");

    % tangent force calculation
    % tangential velocity calculation
    tangvec = [0 -1;1 0]*normvec;
    vecDeltaDot = Rj*(dvecdot/vecnorm(dvec) + ddeltaS*dvec/vecnorm(dvec).^2)-dvecdot;
    vt = dot(vecDeltaDot,tangvec);

    Ft = muf*Fn*min(1,max(0,(abs(vt)-nu0)./(nu1-nu0)))*sign(vt);

    F = Fn*normvec + Ft*tangvec;


end
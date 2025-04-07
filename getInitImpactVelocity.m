function preDeltaDot = getInitImpactVelocity(disBR,betaS,qB,dqB,eBJ,deBJ)
    % function for calculating initial impact velocity \dot{\delta}^{-}
    %
    % Input arguments:
    %   disBR -- discrete radius of bearing contour
    %   betaS -- contact angle of contact point in bearing contour
    %   qB -- generalized coordinates of bearing body
    %   dqB -- generalized velocity of bearing body
    %   eBJ -- eccentricity vector from bearing to journal (GCS)
    %   deBJ -- eccentricity vector dot from bearing to journal (GCS)
    
    arguments (Input)
        disBR (:,1) double {mustBeFinite,mustBeNonNan}
        betaS (1,1) double {mustBePositive,mustBeInRange(betaS,0,6.3)}
        qB (4,1) double {mustBeFinite,mustBeNonNan}
        dqB (4,1) double {mustBeFinite,mustBeNonNan}
        eBJ (2,1) double {mustBeFinite,mustBeNonNan}
        deBJ (2,1) double {mustBeFinite,mustBeNonNan}
    end

    arguments (Output)
        preDeltaDot (1,1) double {mustBeFinite,mustBeNonNan}
    end

    % obtain bearing contour from disBR
    N = numel(disBR);
    betas = 0:2*pi/N:2*pi-2*pi/N;
    Rb = @(x) spline([betas,2*pi],[0;disBR;disBR(1);0],mod(x,2*pi));

    vecDB = consTranMat2D(Rb(betaS).*[cos(betaS);sin(betaS)],"vector")*qB - eBJ;
    vecDBdot = consTranMat2D(Rb(betaS).*[cos(betaS);sin(betaS)],"vector")*dqB - deBJ;

    preDeltaDot = -dot(normalize(vecDB,"norm"),vecDBdot);
end
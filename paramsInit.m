function params = paramsInit()
    % Initial parameter configuration for landing gear retraction mechanism
    % Output arguments:
    %   params--> parameters struct with dynamics and kinematic parameters
    %
    % Remarks: the units of all physical parameters are in SI unit (m-kg-s)

    dyn = struct; % dynamics parameters struct
    dyn.M = [7.66,23.51,147.94,39.16,140.77];
    dyn.J = [2.422e-2,1.137,12.98,4.291,10.76];
    dyn.g = 9.8; % gravity constant [unit]

    kine = struct; % kinematic parameters struct
    kine.H = [1.6,0.4,1.495,0.605,0.3,0.6]; % [1.6 0.4 1.495 0.605 0.3 0.6]
    kine.Ls2 = 0.35; % 0.35
    kine.Ls4 = 0.54; % 0.54
    kine.Ls31 = 0.514; % 0.514
    kine.Ls32 = 0.723; % 0.723
    kine.Ls33 = 0.490; % 0.490
    kine.alpha = deg2rad([16.31,24.48,52.05,15.85,22.55]); % [16.31 24.48 ...
                                                           % 52.05 15.85 22.55]
                                                                            

    % parameters of clearance joint A and joint B

    jointA = struct;
    jointA.ce = .9;
    jointA.muf = .2;
    jointA.nu0 = 1e-3;
    jointA.nu1 = 5e-2;
    jointA.Eb = 200e9;
    jointA.Es = 200e9;
    jointA.nub = 0.3;
    jointA.nus = 0.3;
    jointA.Rb = 10e-3;
    jointA.Rj = 9.95e-3;
    jointA.hB = (1 - jointA.nub^2)/(pi*jointA.Eb);
    jointA.hP = (1 - jointA.nus^2)/(pi*jointA.Es);
    jointA.K = 4/(3*pi*(jointA.hB + jointA.hP))*sqrt(jointA.Rb*jointA.Rj/(jointA.Rb-jointA.Rj));
    jointA.Lc = 15e-3;
    jointA.Kw = 8e-14; % Archard coefficient unit [Pa^-1]
    jointA.N = 3600; % Number of discrete points for bearing contour

    jointB = struct;
    jointB.ce = 0.9;
    jointB.muf = 0.2;
    jointB.nu0 = 1e-3;
    jointB.nu1 = 5e-2;
    jointB.Eb = 200e9;
    jointB.Es = 200e9;
    jointB.nub = 0.3;
    jointB.nus = 0.3;
    jointB.Rb = 10e-3;
    jointB.Rj = 9.95e-3;
    jointB.hB = (1 - jointB.nub^2)/(pi*jointB.Eb);
    jointB.hP = (1 - jointB.nus^2)/(pi*jointB.Es);
    jointB.K = 4/(3*pi*(jointB.hB+jointB.hP))*sqrt(jointB.Rb*jointB.Rj/(jointB.Rb-jointB.Rj));
    jointB.Lc = 15e-3;
    jointB.Kw = 8e-14; % Archard coefficient unit [Pa^-1]
    jointB.N = 3600; % Number of discrete points for bearing contour

    params = struct;
    params.dyn = dyn;
    params.kine = kine;
    params.jointA = jointA;
    params.jointB = jointB;
end
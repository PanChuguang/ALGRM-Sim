function C = consTranMat2D(locPV,type)
    % ----------------------- %
    % constant transform matrix C=[1 0 Px -Py;0 1 Py Px] for point P
    % constant transform matrix C=[0 0 Vx -Vy;0 0 Vy Vx] for vector V
    % for more details, refer to reference ...
    % DOI: 10.1016/j.mechmachtheory.2022.105134
    %
    % Input arguments: locPV --- localized coordinates of point or vector
    %                  type --- Type string
    %
    % Output arguments: C --- transformation matrix from q to gobal,...
    %                         coordinates

    arguments (Input)
        locPV (2,1) {mustBeNumeric}
        type {mustBeMember(type,["point","vector"])}="point"
    end
    
    arguments (Output)
        C (2,4) {mustBeNumeric}
    end

    if type=="point"
        C=[1 0 locPV(1) -locPV(2);0 1 locPV(2) locPV(1)];
    end

    if type=="vector"
        C=[0 0 locPV(1) -locPV(2);0 0 locPV(2) locPV(1)];
    end

end


function [nextq,nextqdot,nextqddot] = oneStepRK4Int(qn,qdn,n,dt,fcnsCell)
    % One time step integration of MBD EoM using RK4 method.
    % 
    % Input arguments:
    %   qn ==> current generalized position
    %   qdn ==> current generalized velocity
    %   n ==> current number of step {tn = (n-1)*dt}
    %   dt ==> fixed time step between adjacent time instants
    %   fcnsCell ==> cell array of fuction handles for DAEs solving
    %     fcnsCell={M,Phi,Phi2q,Phi2qt,Phi2t,Phi2tt,Phi2qqd2q,Q} ordering
    %     signature: M(q),Phi(q,t),Phi2q(q,t),Phi2qt(q,t),Phi2t(q,t),...
    %                Phi2tt(q,t),Phi2qqd2q(q,qdot,t),Q(q,qdot,t)
    %  
    %
    % Output arguments:
    %   nextq <== next time step generalized position 
    %   nextqdot <== next time step generalized velocity
    %
    
    % unpack function handles of DAEs from input argument fcnsCell.
    assert(numel(fcnsCell)==8,"number of function handle must be 8.");
    M = fcnsCell{1};
    Phi = fcnsCell{2};
    Phi2q = fcnsCell{3};
    Phi2qt = fcnsCell{4};
    Phi2t = fcnsCell{5};
    Phi2tt = fcnsCell{6};
    Phi2qqd2q = fcnsCell{7};
    Q = fcnsCell{8};

    % 
    tn = (n-1)*dt;
    numq = length(qn);
    numc = size(Phi(qn,tn),1);

    % RK4 method
    Yn = [qn;qdn];
    k1 = F(Yn,tn);
    k2 = F(Yn+dt/2*k1,tn+dt/2);
    k3 = F(Yn+dt/2*k2,tn+dt/2);
    k4 = F(Yn+dt*k3,tn+dt);
    Yn1 = Yn + dt/6*(k1+2*k2+2*k3+k4);
    nextq = Yn1(1:numq);
    nextqdot = Yn1(numq+1:end);
    Ydotn1 = F(Yn1,tn+dt);
    nextqddot = Ydotn1(numq+1:end);

    
    function Ydot = F(Y,t)
        % ODE for RK4 integration
        % Y = [q;qdot];
        % Ydot = [qdot;qddot]
        
        q = Y(1:numq);
        qdot = Y(numq+1:end);

        Ydot = zeros(2*numq,1);
        Ydot(1:numq) = qdot;

        LHS = [M(q) Phi2q(q,t).';Phi2q(q,t) zeros(numc)];
        RHS = [Q(q,qdot,t);-Phi2qqd2q(q,qdot,t)*qdot-2*Phi2qt(q,t)*qdot-Phi2tt(q,t)-...
            2*20*(Phi2t(q,t)+Phi2q(q,t)*qdot)-20^2*Phi(q,t)];
        accLamda = LHS\RHS;
        Ydot(numq+1:end) = accLamda(1:numq);
    end

end


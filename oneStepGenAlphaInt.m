function [qn1,qdn1,qddn1,lamdan1,an1]=oneStepGenAlphaInt(dt,n,qn,qdn,qddn,an,lamdan,fcnsCell)
    % one time step from n to n+1 using Generalized-\alpha integration method
    %
    % Input arguments:
    %   dt  --> time step between n and n+1
    %   n   --> current number of step
    %   qn  --> current generalized position
    %   qdn --> current generalized velocity
    %   qddn --> current generalized acceleration
    %   an  --> current auxiliary acceleration for generalzied-\alpha
    %   lamdan --> current Lagrange Multiplier
    %   *fcnsCell* --> cell array of function handles for DAEs solving
    %   
    %   fcnsCell={M,Phi,Phi2q,Phi2qt,Phi2t,Phi2tt,Phi2qqd2q,Q} follows the order
    %     signature:
    %     M(q),Phi(q,t),Phi2q(q,t),Phi2qt(q,t),Phi2t(q,t),Phi2tt(q,t),... 
    %     Phi2qqd2q(q,qdot,t),Q(q,qdot,t).
    %   Output arguments:
    %   qn1 --> next time step generalized pos
    %   qdn1 --> next time step generalized vel
    %   qddn1 --> next time step generalized acc
    %   lamdan1 --> next time step Lagrange Multiplier
    %
    
    % get number of generalized coordinates
    ndim = length(qn);
    
    % unpack function handles of DAEs from input arguments *fcnsCell*
    assert(numel(fcnsCell)==8,"The number of function handles must be 8!");
    Mfcn = fcnsCell{1};
    Phifcn = fcnsCell{2};
    Phi2qfcn = fcnsCell{3};
    Phi2qtfcn = fcnsCell{4};
    Phi2tfcn = fcnsCell{5};
    Phi2ttfcn = fcnsCell{6};
    Phi2qqd2qfcn = fcnsCell{7};
    Qfcn = fcnsCell{8};

    %% generalized-\alpha method parameters
    rho = 1;
    am = (2*rho - 1)./(rho + 1);
    af = rho./(rho + 1);
    beta = 1/4*(1 + af - am).^2;
    gam = 1/2 + af - am;

    %% solving nonlinear residue equations r_{n+1} = 0 using "fsolve"
    options = optimoptions("fsolve","Algorithm","levenberg-marquardt","Display","none");
    resfcn = @(X) residueEqns(X(1:ndim),X(ndim+1:end));
    Xn1 = fsolve(resfcn,[an;lamdan],options);
    
    % obtain the next time step quantities
    qn1 = qn + dt*qdn + dt^2*((1/2 - beta)*an + beta*Xn1(1:ndim));
    qdn1 = qdn + dt*(1-gam)*an + gam*dt*Xn1(1:ndim);
    qddn1 = (1-am)/(1-af)*Xn1(1:ndim) + am/(1-af)*an - af/(1-af)*qddn;
    lamdan1 = Xn1(ndim+1:end);
    an1 = Xn1(1:ndim);

    function res=residueEqns(an1,ln1)
        % ------------------ %
        % input arguments:
        %   an1 --> a_{n+1} 
        %   ln1 --> \lamda_{n+1}
        %   
        %  output arguments:
        %   res --> residue equation r_{n+1}
        % 
        
        qnp1 = qn + dt*qdn + dt^2*((1/2 - beta)*an + beta*an1);
        qdnp1 = qdn + dt*(1-gam)*an + gam*dt*an1;
        qddnp1 = (1-am)/(1-af)*an1 + am/(1-af)*an - af/(1-af)*qddn;

        res=[Mfcn(qnp1)*qddnp1+Phi2qfcn(qnp1,n*dt).'*ln1 - Qfcn(qnp1,qdnp1,n*dt);
             Phi2qfcn(qnp1,n*dt)*qddnp1 + Phi2qqd2qfcn(qnp1,qdnp1,n*dt)*qdnp1 + 2*Phi2qtfcn(qnp1,n*dt)*qdnp1 + Phi2ttfcn(qnp1,n*dt)+ ...
             2*20*(Phi2tfcn(qnp1,n*dt)+Phi2qfcn(qnp1,n*dt)*qdnp1)+20^2*Phifcn(qnp1,n*dt)]; % Baumgarte’s stabilization method

    end
   
end
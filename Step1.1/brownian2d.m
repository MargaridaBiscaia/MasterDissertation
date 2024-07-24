function [t,w1,w2] = brownian2d(rho,n)
    % 
    % rho   correlation
    %
    dt = 1/(n-1);
    t = linspace(0,1,n);
    dw = mvnrnd([0.0;0.0],[dt rho*dt;rho*dt dt],n);
    dw(1,:) = [0 0];
    w1 = cumsum(dw(:,1));
    w2 = cumsum(dw(:,2));
    %
end
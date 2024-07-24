function [t,S,V] = hestonSimul(r, S0, varsigma, kappa, delta, v0, rho, n)
    %
    % delta     volatility of volatility
    % rho       correlation
    % kappa     mean reversion speed
    % varsigma  long-term variance
    % v0        initial variance  
    %
    [t,w1,w2] = brownian2d(rho,n+1);
    dt = 1/n;
    V = zeros(n+1,1);
    S = zeros(n+1,1);
    V(1) = v0;
    S(1) = S0;
    %
    for i=2:(n+1)
        dw1 = w1(i) - w1(i-1);
        dS = r*S(i-1)*dt + sqrt(v0)*S(i-1)*dw1;
        S(i) = S(i-1) + dS;
        dw2 = w2(i) - w2(i-1);
        dV = kappa*(varsigma - v0)*dt + delta*sqrt(v0)*dw2;
        v0 = v0 + dV;
        V(i) = v0;
    end
    %
    t = t';
    S = round(real(S),2);
    %
    %
end
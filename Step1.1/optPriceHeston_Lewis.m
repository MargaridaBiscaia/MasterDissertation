function [price] = optPriceHeston_Lewis(K, S, r, tau, delta, rho, kappa, varsigma, v0)
    %
    % delta     volatility of volatility
    % rho       correlation
    % kappa     mean reversion speed
    % varsigma  long-term variance
    % v0        initial variance  
    %
    function cf = CF(u, j, S, r, tau, delta, rho, kappa, varsigma, v0)
        %
        a = [0.5, -0.5];
        b = [kappa - rho*delta, kappa];
        %
        d = sqrt((b(j) - rho*delta*1i*u).^2 - delta^2*(2*a(j)*1i*u - u.^2));
        g = (b(j) - rho*delta*1i*u-d) / (b(j) - rho*delta*1i*u+d);
        %
        C = 1i*u*r*tau + ((varsigma*kappa)/delta^2)*((b(j) - rho*delta*1i*u -d)*tau - 2*log((1-g*exp(-d*tau))/(1-g)));
        D = (b(j) - rho*delta*1i*u - d)/(delta^2) * (1-exp(-d*tau))/(1-g*exp(-d*tau));
        %
        cf = exp(C + D*v0 + 1i*u*log(S));
        %
    end
    %
    cf2 = @(u) CF(u, 2, S, r, tau, delta, rho, kappa, varsigma, v0);
    i = @(u)  K^(-1i*u) * cf2(u - 1i/2) / (u^2 + 0.25) ;
    ii = @(u) real(i(u));
    
    tt =  integral(@(u) ii(u), 0, inf, 'ArrayValued', true);
    %
    price = S - sqrt(K)*exp(-r * tau) / pi * tt;
    %
end

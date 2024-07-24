%% Heston %% 
Np = 1; % number of sets of parameters
D = 5; % number of parameters
%
n = 52; % number of observations (1 observation per day in a year)
r = 0.05; % risk-free interest rate
S0 = 189; % initial stock price (thinking of apple)
%
Nop = 550; % number of options
%
%% Generate Heston parameters
varsigma = 0.1;
kappa = 1;
delta = 0.2;
v0 = 0.04;
rho = -0.75;
parameters = [varsigma, kappa, delta, v0, rho];
%% Calculate correspondent stock prices, for each set of parameters
matrizStock = zeros(n+1, Np);
%
for i=1:Np
    [t, S] = hestonSimul(r, S0, varsigma(i), kappa(i), delta(i), v0(i), rho(i), n);
    while any(S < 164) || any(S > 198)
        [t, S] = hestonSimul(r, S0, varsigma(i), kappa(i), delta(i), v0(i), rho(i), n);
    end

    matrizStock(:,i) = S(:);
    %
    plot(t, S);
    xlabel('Time, t (years)');
    ylabel('Stock Price, St');
    title('Stock Prices with the Heston Model');
    hold on;
end
hold off;
%
%% Replicate the parameters (n+1) times (goal: parameters organized per day) and associate the stock prices to each line
dataset1 = [repmat(parameters, n+1, 1), reshape(matrizStock.', [], 1)];
%
%% Generate derivatives properties (Maturity time T and Strike prices)
T = linspace(1, 2.0, 11)';
Ns = Nop/11; % Number of strikes per maturity
strike = linspace(151, 249, Ns)';
%%
%
options = [repelem(T, Ns,1), repmat(sort(strike), 11, 1)];
%
%% Associate each line to an option 
dataset2 = [repelem(dataset1,Nop,1), repmat(options,(n+1)*Np,1)];
%
%% Determine tau - time to maturity (if the data is daily!!!!)
dataset3 = dataset2;
%
for i=1:(Nop*Np)
    dataset3(i,9) = dataset2(i,7);
end
%
for k=1:n
    for i = (Nop*Np*k+1):(Nop*Np*k+Nop*Np)
        dataset3(i,9) = abs(dataset2(i,7)-k/52);
    end
end
%
dataset3 = array2table(dataset3);
dataset3.Properties.VariableNames = {'varsigma', 'kappa', 'delta', 'v0', 'rho', 'Stock Price', 'Maturity Time', 'Strike', 'tau'};
%
paraAqui
%% Calculate the prices
ntotal = Nop*Np*(n+1);
Price = zeros(ntotal, 1);
date = 0;
Maturity = floor(dataset3.tau*365+1);
%%
%
h = waitbar(0,'Please wait ...');
niter = ntotal;
iter = 0;

for i = 1:ntotal
    pricec = optByHestonNI(r, dataset3.('Stock Price')(i), date, Maturity(i), 'call', dataset3.Strike(i), dataset3.v0(i), dataset3.varsigma(i), dataset3.kappa(i), dataset3.delta(i), dataset3.rho(i), 'DividendYield', 0, 'Framework', 'lewis2001');
    iter = iter + 1;
    waitbar(iter/niter,h,'Please wait ...')
    %pricec = optPriceHeston_Lewis(dataset3.Strike(i), dataset3.('Stock Price')(i), r, dataset3.tau(i), dataset3.delta(i), dataset3.rho(i), dataset3.kappa(i), dataset3.varsigma(i), dataset3.v0(i));
    Price(i) = pricec;
end 
close(h)
%
j=i;
save('Price.mat', 'Price', 'j');
%
% Calculate moneyness and obtain the final dataset
moneyness = dataset3.('Stock Price') ./ dataset3.Strike;
%
% Final dataset
dataset = [dataset3(:, {'varsigma', 'kappa', 'delta', 'v0', 'rho', 'tau', 'Stock Price', 'Strike'}), array2table(moneyness), array2table(Price)];
dataset = round(table2array(dataset),4);
dataset = array2table(dataset);
dataset.Properties.VariableNames = {'varsigma', 'kappa', 'delta', 'v0', 'rho', 'tau', 'stockPrice', 'strike', 'moneyness', 'price'};
%
writetable(dataset,"C:\Users\User\OneDrive - Universidade de Coimbra\Ambiente de Trabalho\Thesis\Code\Apple\calibrateApple.csv")
writetable(dataset3,"C:\Users\User\OneDrive - Universidade de Coimbra\Ambiente de Trabalho\Thesis\Code\Apple\featuresCalibrateApple.csv")
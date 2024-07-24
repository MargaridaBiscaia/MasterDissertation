%% Heston %% 
Np = 300; % number of sets of parameters
D = 5; % number of parameters
%
n = 52; % number of observations (1 observation per week in a year)
r = 0.05; % risk-free interest rate
S0 = 180; % initial stock price (thinking of apple)
%
Nop = 110; % number of options
%
%% Generate Heston parameters
% Define the ranges for each parameter
range_varsigma = [0.01, 0.5]; % long-term variance
range_kappa = [0, 3.0];
range_v0 = [0.03, 0.15];
range_delta = [0.01, 0.8]; % volatility of volatility
range_rho = [-0.9, -0.2];
%% Generate Latin Hypercube Samples
X = lhsdesign(Np, D, 'criterion', 'maximin', 'iterations', 500);

% Rescale the samples to the specified ranges
varsigma = range_varsigma(1) + X(:, 1) * (range_varsigma(2) - range_varsigma(1));
kappa = range_kappa(1) + X(:, 2) * (range_kappa(2) - range_kappa(1));
v0 = range_v0(1) + X(:, 3) * (range_v0(2) - range_v0(1));
delta = range_delta(1) + X(:, 4) * (range_delta(2) - range_delta(1));
rho = range_rho(1) + X(:, 5) * (range_rho(2) - range_rho(1));

% Ensure the condition 2*kappa*varsigma > delta^2 is met
valid_indices = 2 * kappa .* varsigma > delta .^ 2;
while ~all(valid_indices)
    % Regenerate the samples for the invalid indices
    invalid_X = lhsdesign(sum(~valid_indices), D, 'criterion', 'maximin', 'iterations', 500);
    varsigma(~valid_indices) = range_varsigma(1) + invalid_X(:, 1) * (range_varsigma(2) - range_varsigma(1));
    kappa(~valid_indices) = range_kappa(1) + invalid_X(:, 2) * (range_kappa(2) - range_kappa(1));
    v0(~valid_indices) = range_v0(1) + invalid_X(:, 3) * (range_v0(2) - range_v0(1));
    delta(~valid_indices) = range_delta(1) + invalid_X(:, 4) * (range_delta(2) - range_delta(1));
    rho(~valid_indices) = range_rho(1) + invalid_X(:, 5) * (range_rho(2) - range_rho(1));
    
    % Revalidate the condition
    valid_indices = 2 * kappa .* varsigma > delta .^ 2;
end

% Round the parameters to 2 decimal places
varsigma = round(varsigma, 2);
kappa = round(kappa, 2);
v0 = round(v0, 2);
delta = round(delta, 2);
rho = round(rho, 2);

% Create the parameters matrix
parameters = [varsigma, kappa, delta, v0, rho];
%
%% Calculate correspondent stock prices, for each set of parameters
matrizStock = zeros(n+1, Np);
%
for i=1:Np
    [t, S] = hestonSimul(r, S0, varsigma(i), kappa(i), delta(i), v0(i), rho(i), n);
    while any(S < 150) || any(S > 210)
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
strike = [150, 160, 170, 180, 190, 200, 210, 220, 230, 250]';
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
%% Calculate the prices
ntotal = Nop*Np*(n+1);
Price = zeros(ntotal, 1);
%%
%
h = waitbar(0,'Please wait ...');
niter = ntotal;
iter = 0;

for i = 1:ntotal
    iter = iter + 1;
    waitbar(iter/niter,h,'Please wait ...')
    pricec = optPriceHeston_Lewis(dataset3.Strike(i), dataset3.('Stock Price')(i), r, dataset3.tau(i), dataset3.delta(i), dataset3.rho(i), dataset3.kappa(i), dataset3.varsigma(i), dataset3.v0(i));
    Price(i) = pricec;
end 
close(h)
%
% Calculate moneyness and obtain the final dataset
moneyness = dataset3.('Stock Price') ./ dataset3.Strike;
%
% Final dataset
dataset = [dataset3(:, {'varsigma', 'kappa', 'delta', 'v0', 'rho', 'tau', 'Stock Price', 'Strike'}), array2table(moneyness), array2table(Price)];
dataset = round(table2array(dataset),4);
dataset = array2table(dataset);
dataset.Properties.VariableNames = {'varsigma', 'kappa', 'delta', 'v0', 'rho', 'tau', 'stockPrice', 'strike', 'moneyness', 'price'};

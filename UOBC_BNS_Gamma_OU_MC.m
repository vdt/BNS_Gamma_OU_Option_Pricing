% Ioannis Anagnostou, 2010

% function UOBC_BNS_Gamma_OU_MC uses Monte Carlo simulation to
% estimate the price of a Up-and-out Barrier Call option under the BNS model
% with gamma-OU stochastic volatility

% S0: initial price of the underlying
% K: option strike
% H: barrier
% r: annualized risk-free interest rate
% q: divident
% sigma: volatility of the underlying
% T: time to maturity
% N: time steps
% n: number of simulations
% lambda, a, b, and rho: parameters of the BNS Gamma-OU model

function UOBC_BNS_Gamma_OU_MC(lambda,a,b,rho,S0,K,H,r,q,sigma,T,N,n)

dt = T/N; t = (0:dt:T);
I = zeros(n,N);

PS = zeros(n,N+1);
Y = zeros(n,N+1);
Z = zeros(n,N+1);
BNS = zeros(n,N+1);

for i = 1:n

    for j = 1:N
        I(i,j) = poissrnd(dt*a*lambda);
        PS(i,j+1) = PS(i,j) + I(i,j);
    end

    for j = 2:N
        if PS(i,j) > PS(i,j-1)
            dif = PS(i,j)-PS(i,j-1);
            S = zeros(dif+1,1);
            for k = 1:dif
            S(k) = Gamma(1,b);
            end
            C = cumsum(S);
            Y(i,j) = C(dif,:);
        else
            Y(i,j) = 0;
        end
    end

    Z(i,1) = 0.01452;
    
    for j = 1:N
        Z(i,j+1) = (1-lambda*dt)*Z(i,j)+Y(i,j+1);
    end

    for j = 1:N
        BNS(i,j+1) = BNS(i,j)+(r-q-((lambda*rho*a)/(b-rho))-(Z(i,j+1)/2))*dt+sqrt(Z(i,j+1)*dt)*randn+rho*(Y(i,j+1)-Y(i,j));
    end

end

S = S0*exp(BNS);
ST = S(:,N+1);
m = max(S,[],2); 

for i = 1:n  
    if m(i) <= H
        c(i) = max(0,ST(i)-K);
    else 
        c(i) = 0;
    end
end    
    
call_price = exp(-r*T)*mean(c)
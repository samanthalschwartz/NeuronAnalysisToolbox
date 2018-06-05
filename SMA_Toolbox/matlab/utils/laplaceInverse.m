%Mark J. Olah (mjo@cs.unm.edu)
% 06-2016
% Numerical laplace inverse with the Gaver-Stehfest algorithm

function val = laplaceInverse(F, t, M)
    % [in]
    %  F - a function handle that accepts 's' laplace space coordinates
    %  t - time to evaluate at
    %  M - [optional] Order (default=7)
    
    if nargin<3
        M=7;
    end
    Q = Qfactors(M);
    Qsum = sum(Q,2);
    omega = zeros(2*M,1);
    Fvals = zeros(2*M,1);

    alpha = log(2)/t;
    for k=1:2*M
        omega(k) = (-1)^(M+k) * Qsum(k);
        Fvals(k) = omega(k)*F(k*alpha);
    end
    val = alpha * sum(Fvals);
end

function Q = Qfactors(M)
    Q=zeros(2*M,M);
    for k=1:2*M
        for j = floor((k+1)/2):min(k,M)
            Q(k,j) = j^(M+1) / (factorial(j)^2*factorial(M-j)*factorial(k-j)*factorial(2*j-k)) * factorial(2*j);
        end
    end
end
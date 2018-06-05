% Mark J. Olah [mjo@cs.unm.edu]
% 2016

function norm_prob = logNormalize(log_likelihood,MIN_LOG)
    % [in]
    %  log_likelihood - size:[N] array of log-liklihoods
    %  MIN_LOG - [optional] Minimum log value to exponenitate (smaller values are 0)  [default: -1e2]
    % [out]
    %  norm_prob - size:[N] array of normalized probabilities
    %  logSum - log of normalization constant. Equivalent to: log(sum(exp(log_likelihood))) but computed as accurately as possible    
    if nargin<2
        MIN_LOG = -1e2;
    end
    log_likelihood = log_likelihood - max(log_likelihood);
    norm_prob = zeros(size(log_likelihood));
    non_zero = log_likelihood>MIN_LOG;
    norm_prob(non_zero) = exp(log_likelihood(non_zero));
    norm_prob = norm_prob/sum(norm_prob);
end

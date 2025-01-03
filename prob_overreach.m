function [prob_error_sum ] = prob_overreach(distances, noise_ref_linear)
%given intersymbol distances, computes the probability of symbol misdetection
    N = numel(distances);
    prob_nothing_goes_wrong = 1;
    for k = 2:1:N
        d = distances(k);
        prob_k = qfunc((d*0.5)/((0.5*noise_ref_linear)^0.5));
        inv_prob = 1 - prob_k;
        prob_nothing_goes_wrong = prob_nothing_goes_wrong * inv_prob;
        
    end
    prob_something_goes_wrong = 1 - prob_nothing_goes_wrong;
    prob_error_sum = prob_something_goes_wrong;

end
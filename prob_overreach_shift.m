function [prob_error_sum, prob_good_sum] = prob_overreach_shift(distance_mtx, noise_ref_linear)
    N = size(distance_mtx, 1);
    %NOTE: the distance mtx is a 2D matrix. Diagonal elements contain the
    %distance between the original elements and their shifted variants,
    %while non-diag elements are distances to other point in constellaion
    %diagram
    prob_dissimilar_mismatch = 1;
    prob_similar_mismatch = 1;
    
    for i = 1:1:N
        for j = 1:1:N
            d = distance_mtx(i,j);
            prob_k = qfunc((d*0.5)/((0.5*noise_ref_linear)^0.5));
            if (i ~= j)
                prob_dissimilar_mismatch = prob_dissimilar_mismatch * (1 - prob_k); % aggregate probability of symbols being mistaken for ezach other
            else %two same symbols
                prob_similar_mismatch = prob_similar_mismatch * (1 - prob_k);

            end
        end
    end
    %note: symbol is misrecognized if it is misttributed to other symbol OR
    %not associated with its non shifted self
        prob_error_sum_ex = 1 - prob_dissimilar_mismatch;
        prob_good_sum = 1 - prob_similar_mismatch; 
        remainder = 1 - prob_error_sum_ex - prob_good_sum;
        prob_error_sum = prob_error_sum_ex;% + remainder; %fixme later
    
end
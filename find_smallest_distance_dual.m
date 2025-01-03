function [distance_matrix, nearest_neightbours_different] = find_smallest_distance_dual(pointlist_ref, pointlist_shifted)
%given the constellation and its hsifted variant, computes the distance to
%each similar and non similar point and returns as a matrix
    N1 = numel(pointlist_ref);

    closest_idx_non_analogous = zeros(N1,1);
    dists = zeros(N1,N1);
    for k = 1:1:N1 
        win_idx1 = 0;
        win_idx2 = 0;
        win_dist = inf;
        for j = 1:1:N1 
            p1 = pointlist_ref(k);
            p2 = pointlist_shifted(j);
            p1_r = real(p1);
            p2_r = real(p2);
            p1_i = imag(p1);
            p2_i = imag(p2);
            distance_21 = (p2_i - p1_i)^2 + (p2_r - p1_r)^2;
            distance_21 = distance_21 ^0.5;
            %if k == j analogous points, matched symbols
            dists(k,j) = distance_21;
            if (win_dist > distance_21) && (k ~= j)
                win_dist = distance_21;
                closest_idx_non_analogous(k) = j;
            end

        end
    end
    distance_matrix = dists;
    nearest_neightbours_different = closest_idx_non_analogous;

end
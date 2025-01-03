function [d, idx1, idx2] = find_smallest_distance_complex(pointlist)
    %returns, given an array of imaginary/2d numbers, for each element, the
    % index of the closest point and the distance
    N = numel(pointlist);
    idx1_l = zeros(N, 1);
    idx2_l = zeros(N,1);
    win_dists = zeros(N,1);


    for u = 1:1:N-1 
        win_idx1 = 0;
        win_idx2 = 0;
        win_dist = inf;
        
        for v = u+1:1:N 
            p1 = pointlist(u);
            p2 = pointlist(v);
            p1_r = real(p1);
            p2_r = real(p2);
            p1_i = imag(p1);
            p2_i = imag(p2);

            %euclidian distance between two symbols
            distance_21 = (p2_i - p1_i)^2 + (p2_r - p1_r)^2;
            distance_21 = distance_21 ^0.5;
            if distance_21 < win_dist %if new winner found
                win_idx1 = v;
                win_idx2 = u;
                win_dist = distance_21;
            
            end
        idx1_l(u) = win_idx1;
        idx2_l(u) = win_idx2;
        win_dists(u) = win_dist;

        end
    end
    P = [win_dists, idx1_l, idx2_l];
    P = sortrows(P);

    d = P(:,1);
    idx1 = P(:,2);
    idx2 = P(:,3);
end
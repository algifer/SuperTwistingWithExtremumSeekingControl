function lhs_points = generate_lhs_points(n_samples, n_dimensions, bounds)
    % Genera puntos utilizando Latin Hypercube Sampling (LHS)
    lhs_points = lhsdesign(n_samples, n_dimensions);
    scaled_points = zeros(size(lhs_points));
    
    for i = 1:n_dimensions
        lower = bounds(i, 1);
        upper = bounds(i, 2);
        scaled_points(:, i) = lower + (upper - lower) * lhs_points(:, i);
    end
    
    lhs_points = scaled_points;
end
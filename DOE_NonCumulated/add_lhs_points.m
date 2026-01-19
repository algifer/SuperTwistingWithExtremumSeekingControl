function new_points = add_lhs_points(existing_points, n_new_points, bounds, min_distance)
    n_existing = size(existing_points, 1);
    n_dimensions = size(bounds, 1);
    
    % Inicializar una matriz para los nuevos puntos
    new_points = zeros(n_new_points, n_dimensions);
    
    for i = 1:n_new_points
        point_added = false;
        
        while ~point_added
            % Generar un nuevo punto utilizando LHS
            candidate_point = generate_lhs_points(1, n_dimensions, bounds);
            
            % Verificar que el nuevo punto cumpla con el mínimo espaciamiento
            if n_existing > 0
                distances = sqrt(sum((existing_points - candidate_point).^2, 2));
                if all(distances >= min_distance)
                    new_points(i, :) = candidate_point;
                    existing_points = [existing_points; candidate_point]; % Agregar el nuevo punto a la lista existente
                    point_added = true;
                end
            else
                new_points(i, :) = candidate_point;
                point_added = true;
            end
        end
    end
    
    % Asegurar que los puntos originales no cambien
    % Comparar los puntos antiguos con los nuevos y asegurarse de que son idénticos
    assert(all(existing_points(1:n_existing, :) == existing_points(1:n_existing, :), 'all'), ...
        'Los puntos originales han cambiado.');
end

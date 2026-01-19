%DOE Puntos Nuevos
function main()
    % Definir los límites para cada dimensión (en este caso, dos dimensiones)
    bounds = [0.1 40; 0.1 5.5];
    end_points = 30;
    
    initial_points = [
        16.7772,    3.8975;
       14.9453,    0.8458;
       36.1393,    3.1449;
       28.3767,    4.5859;
        2.1762,    1.9429;
       12.1258,    1.6992;
       27.7237,    3.9473;
        1.7015,    5.4231;
       36.7284,    0.8911;
       17.1324,    3.0863;
       21.5217,    2.4688;
        5.7501,    5.3628;
       38.1810,    4.3480;
       10.2701,    2.0390;
       29.1764,    1.0099;
        1.3451,    1.0355;
       26.4169,    5.1796;
       15.8845,    1.4789;
       18.7638,    3.4218;
       38.9407,    3.2193;
       34.9423,    2.2638;
        3.5327,    1.4701;
       22.7139,    5.1013;
        8.2167,    4.1480;
       31.3988,    0.1868;
        1.5889,    1.8267;
       31.4892,    3.4130;
       38.1897,    0.2909;
       16.3567,    2.3405;
       15.2324,    5.3212;
    ];

    % Verificar si se deben generar nuevos puntos iniciales
    if isempty(initial_points)
        generate_initial_points = true;
        add_new_points = false;
    else
        generate_initial_points = false;
        add_new_points = true;
    end

    % Generar puntos iniciales si está configurado para hacerlo
    if generate_initial_points
        initial_points = generate_lhs_points(5, 2, bounds);
        disp('Puntos iniciales generados:');
        disp(initial_points);  
    else
        disp('Puntos iniciales previos:');
        disp(initial_points);   
    end
    
    if add_new_points && (length(initial_points) < end_points)
        % Agregar nuevos puntos asegurando un mínimo espaciamiento
        new_points = add_lhs_points(initial_points, 5, bounds, 0.2);
        disp('Nuevos puntos agregados:');
        disp(new_points);

        % Combinar puntos iniciales y nuevos puntos
        all_points = [initial_points; new_points];
        disp('Total de puntos:');
        disp(all_points);
    else
        all_points = initial_points;
        disp('Total de puntos:');
        disp(all_points);
    end
    
    if length(initial_points) <= end_points
        % Evaluar puntos y guardar resultados en Excel
        evaluate_and_save_points(all_points, bounds);
    else
        
    end
end

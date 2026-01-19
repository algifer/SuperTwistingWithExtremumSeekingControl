function main()
    clc, close all, clear all; 
    % Definir los límites para cada dimensión (en este caso, dos dimensiones)
    bounds = [0.1 30; 0.1 5.6];
    end_points = 35;
    
    initial_points = [
        6.4777,    4.9908;
        2.1489,    1.7942;
       15.8655,    3.1388;
       28.1051,    3.5445;
       19.6798,    0.3754;
       3.5032,    4.1413;
        1.3241,    2.9175;
       24.9560,    5.2023;
       23.7996,    3.4776;
       23.4657,    0.1005;
       25.8004,    1.7185;
        9.7243,    4.5825;
       11.6560,    0.1604;
       13.3452,    1.6807;
       11.0930,    2.5240;
       12.5836,    5.3628;
       20.3222,    2.1217;
       16.7966,    3.5058;
       13.2007,    3.6085;
       28.3051,    2.6594;
       13.6396,    3.3582;
       28.0157,    2.2984;
       17.5394,    1.1139;
       25.5518,    1.6059;
       18.9806,    1.2686;
       15.3872,    5.2852;
        2.4775,    5.4531;
        4.5742,    4.6849;
       26.6401,    2.1402;
       14.3322,    1.2821;
       17.4552,    4.9028;
       26.5202,    2.3625;
       22.7237,    5.0014;
        0.9306,    1.2633;
        7.7417,    1.3076;
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
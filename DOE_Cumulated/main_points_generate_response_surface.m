function main_points()
    % Definir los límites para cada dimensión (en este caso, dos dimensiones)
    bounds = [0.1 30; 0.1 5.6];
    end_points = 35;
    
    initial_points = [
    4.2621,    2.5632;
    7.6858,    1.4270;
   10.9864,    1.5621;
   19.4166,    3.3720;
   17.9516,    2.9327;
    1.0477,    5.3572;
   14.9297,    1.6943;
   22.9111,    2.7454;
   28.1463,    2.4312;
    5.7695,    2.2866;
    1.9919,    2.8615;
   18.3028,    1.9503;
   20.1162,    5.0819;
    3.4913,    0.8960;
    6.3696,    0.3271;
   16.7475,    2.1096;
    5.0948,    4.7994;
   29.2612,    0.5148;
    7.9837,    5.1720;
   26.0708,    3.4576;
   12.9821,    3.1197;
   23.7329,    4.2691;
   26.6920,    1.1714;
   20.9033,    3.6392;
   25.4083,    0.1980;
   28.7738,    4.4555;
    8.8834,    5.5037;
   21.6697,    3.8774;
   24.6279,    1.3459;
   15.7055,    0.8803;
   14.2965,    4.5766;
   12.0965,    3.7814;
   11.2302,    4.0318;
    9.8581,    0.6170;
    0.9275,    4.9143;
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
        initial_points = generate_lhs_points(end_points, 2, bounds);
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
        evaluate_and_save_points_main_points(all_points, bounds);
    else
        
    end
end

%DOE Acumulative
function main()
    clc, close all, clear all; 
    % Definir los límites para cada dimensión (en este caso, dos dimensiones)
    bounds = [0.1 30; 0.1 5.6];
    end_points = 35;
    
    initial_points = [
        24.8386,    3.0337;
       22.0812,    5.5171;
        4.2528,    4.0629;
        7.2261,    0.4728;
       15.7592,    1.3794;
        1.0560,    1.2809;
        8.6776,    1.3766;
        4.0497,    1.2964;
       11.0605,    2.5413;
       18.0539,    3.6159;
       10.5595,    5.0666;
       29.7624,    0.3038;
        5.1193,    4.7701;
        3.0425,    4.8433;
       16.9176,    4.1118;
        5.7014,    1.3615;
       19.7977,    2.0498;
       25.1342,    2.4285;
        5.9654,    1.2211;
       20.5506,    4.0067;
       22.4360,    0.3670;
       19.4721,    3.0439;
       14.7590,    3.7973;
       29.7484,    3.1410;
       24.3058,    3.2698;
       19.9761,    1.1173;
       19.9744,    3.5834;
       18.2835,    3.3613;
       29.2295,    2.1980;
       29.6155,    4.0041;
       28.8648,    1.0760;
       22.1499,    2.7267;
       27.9251,    5.4577;
       19.4611,    3.7806;
       27.9258,    0.2461;
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
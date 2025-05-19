% ======== FUNCIÓN PRINCIPAL PARA PROCESAR GEOMETRÍAS ========
function geo_data = procesar_geometrias()
    % Buscar automáticamente el archivo de salida
    posibles_nombres = {'*.nwo', '*.nwout', '*.out','*.log'};
    archivo_encontrado = '';
    
    for i = 1:length(posibles_nombres)
        archivos = dir(posibles_nombres{i});
        if ~isempty(archivos)
            archivo_encontrado = archivos(1).name;
            break;
        end
    end
    
    if isempty(archivo_encontrado)
        error('No se encontró archivo de salida');
    end
    
    % Procesar con el nombre correcto
    geo_data = procesar_geometrias(archivo_encontrado);
end

% ======== FUNCIONES AUXILIARES ========
function distancias = calcular_distancias(coordenadas)
    n_atoms = size(coordenadas,1);
    distancias = zeros(n_atoms);
    for i = 1:n_atoms
        for j = i+1:n_atoms
            distancias(i,j) = norm(coordenadas(i,:) - coordenadas(j,:));
        end
    end
end

function angulos = calcular_angulos(coordenadas)
    n_atoms = size(coordenadas,1);
    angulos = [];
    for i = 1:n_atoms
        for j = i+1:n_atoms
            for k = j+1:n_atoms
                v1 = coordenadas(j,:) - coordenadas(i,:);
                v2 = coordenadas(k,:) - coordenadas(i,:);
                cos_theta = dot(v1,v2)/(norm(v1)*norm(v2));
                angulos = [angulos; acosd(cos_theta)];
            end
        end
    end
end

function I = calcular_momento_inercia(geo)
    % Masas atómicas (ejemplo para H2O)
    masas = obtener_masas(geo.atomos); 
    
    % Centro de masa
    cm = sum(masas'.*geo.coordenadas)/sum(masas);
    
    % Matriz de inercia
    I = zeros(3);
    for i = 1:length(masas)
        r = geo.coordenadas(i,:) - cm;
        I = I + masas(i)*(r'*r*eye(3) - r'*r);
    end
end

function visualizar_geometria(geo)
    figure;
    scatter3(geo.coordenadas(:,1), geo.coordenadas(:,2), geo.coordenadas(:,3), ...
             200, 'filled');
    text(geo.coordenadas(:,1), geo.coordenadas(:,2), geo.coordenadas(:,3), ...
         geo.atomos, 'FontSize', 14);
    axis equal;
    xlabel('X (Å)'); ylabel('Y (Å)'); zlabel('Z (Å)');
    title('Geometría Molecular');
end

% ======== EJEMPLO DE USO ========
% Cargar y procesar el archivo
geo_data = procesar_geometrias('tu_archivo.nw');

% Acceder a los resultados para la última geometría
ultima_geo = geo_data{end};
disp('Distancias interatómicas (Å):');
disp(ultima_geo.distancias);
disp('Ángulos de enlace (grados):');
disp(ultima_geo.angulos);
disp('Matriz de momento de inercia (u·Å²):');
disp(ultima_geo.momento_inercia);

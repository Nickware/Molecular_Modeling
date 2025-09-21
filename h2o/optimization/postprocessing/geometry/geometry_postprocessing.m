function geo_data = geometry_postprocessing(nombre_archivo)
    % Procesa y analiza geometrías moleculares desde un archivo de salida de NWChem
    geometrias = procesar_geometrias(nombre_archivo);
    for i = 1:length(geometrias)
        geo = geometrias{i};
        geo.distancias = calcular_distancias(geo.coordenadas);
        geo.angulos = calcular_angulos(geo.coordenadas);
        geo.momento_inercia = calcular_momento_inercia(geo);
        geometrias{i} = geo;
    end
    geo_data = geometrias;
    visualizar_geometria(geometrias{end});
end

function geo_data = procesar_geometrias(nombre_archivo)
    % Extrae bloques de geometría molecular del archivo de salida
    geo_data = procesar_geometrias_impl(nombre_archivo);
end

function geo_data = procesar_geometrias_impl(nombre_archivo)
    fid = fopen(nombre_archivo, 'r');
    if fid == -1
        error('Archivo no encontrado');
    end

    geo_data = {};
    leyendo = false;
    coords = [];
    atomos = {};

    while ~feof(fid)
        linea = fgetl(fid);
        % Detectar inicio de bloque de geometría
        if ~isempty(strfind(linea, 'geometry')) && ~isempty(strfind(linea, 'angstroms'))
            leyendo = true;
            coords = [];
            atomos = {};
            continue;
        end

        % Detectar fin de bloque
        if leyendo && (isempty(linea) || ~isempty(strfind(linea, 'end')))
            leyendo = false;
            if ~isempty(coords)
                geo_data{end+1} = struct('atomos', {atomos}, 'coordenadas', coords);
            end
            continue;
        end

        % Procesar líneas de coordenadas
        if leyendo
            tokens = strsplit(strtrim(linea));
            if length(tokens) >= 4
                atomos{end+1} = tokens{1};
                coords = [coords; str2double(tokens(2:4))];
            end
        end
    end

    fclose(fid);
end

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

function masas = obtener_masas(atomos)
    % Tabla simple de masas atómicas para algunos elementos comunes
    tabla_masas = containers.Map({'H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl'}, ...
                                [1.0079, 12.0107, 14.0067, 15.999, 18.998, 30.973, 32.065, 35.453]);
    n = length(atomos);
    masas = zeros(1,n);
    for i = 1:n
        if isKey(tabla_masas, atomos{i})
            masas(i) = tabla_masas(atomos{i});
        else
            masas(i) = 0; % Elemento desconocido
        end
    end
end

function I = calcular_momento_inercia(geo)
    masas = obtener_masas(geo.atomos);
    cm = sum(bsxfun(@times, geo.coordenadas, masas'), 1) / sum(masas);
    I = zeros(3);
    for i = 1:length(masas)
        r = geo.coordenadas(i,:) - cm;
        I = I + masas(i)*(dot(r,r)*eye(3) - (r'*r));
    end
end

function visualizar_geometria(geo)
    figure;
    scatter3(geo.coordenadas(:,1), geo.coordenadas(:,2), geo.coordenadas(:,3), 200, 'filled');
    text(geo.coordenadas(:,1), geo.coordenadas(:,2), geo.coordenadas(:,3), geo.atomos, 'FontSize', 14);
    axis equal;
    xlabel('X (Å)'); ylabel('Y (Å)'); zlabel('Z (Å)');
    title('Geometría Molecular');
end
% parsear_nwchem.m
% Versión robusta para extraer datos de salida de NWChem
% Compatible con diferentes formatos y versiones
% Autor: ChatGPT
% Este script lee un archivo de salida de NWChem, extrae energías DFT, coordenadas atómicas y gradientes,
% y luego calcula distancias y ángulos para la molécula de agua.
% El script es flexible y maneja variaciones en el formato de salida, con mensajes de depuración para facilitar la identificación de problemas.
% Requiere que el archivo "energy.nwo" esté en el mismo directorio.
% Uso:
% 1. Coloca el archivo "energy.nwo" en el mismo directorio que
%    este script.
% 2. Ejecuta el script en Octave o MATLAB.
% El script mostrará un resumen de los datos extraídos y guardará los resultados en "datos_nwchem_extraidos.mat".
% Nota: Asegúrate de tener permisos de lectura para el archivo "energy.nwo" y de que el formato del archivo sea compatible con las expectativas del script.
% Si el script no encuentra energías, coordenadas o gradientes, mostrará advertencias para ayudarte a identificar posibles problemas con el formato del archivo de salida de NWChem.
% El script también incluye gráficos para visualizar la convergencia de la energía, la norma del gradiente, las longitudes de enlace y el ángulo de enlace a lo largo de las iteraciones.
% Si tienes problemas con el formato del archivo, revisa los mensajes de depuración para identificar qué parte del proceso de extracción no está funcionando correctamente.
% Si el formato del archivo de salida de NWChem es muy diferente, es posible que necesites ajustar los patrones de búsqueda en el script para adaptarlo a tu caso específico.
% Este script es un punto de partida robusto para analizar la salida de NWChem, pero puede requerir ajustes dependiendo de la versión de NWChem y del formato específico de tu archivo de salida.
% Si tienes problemas o necesitas ayuda para adaptar el script a tu formato específico, no dudes en pedir asistencia.
% Nota: Este script asume que el archivo de salida de NWChem contiene líneas con la energía total DFT en el formato "Total DFT energy = -76.xxxx Ha", coordenadas atómicas en un bloque que sigue a la línea "Output coordinates in angstroms", y gradientes en un bloque que sigue a la línea "DFT ENERGY GRADIENTS". Si tu archivo de salida tiene un formato diferente, es posible que necesites ajustar los patrones de búsqueda en el script para adaptarlo a tu caso específico.
% ¡Buena suerte con tu análisis de NWChem!
clear; clc;

% Nombre del archivo
filename = "energy.nwo";

if ~exist(filename, "file")
    error("No se encuentra el archivo: %s", filename);
endif

% Leer todo el archivo
fid = fopen(filename, "r");
if fid == -1
    error("No se pudo abrir el archivo: %s", filename);
endif
lines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
lines = lines{1};  % Vector de cadenas

% Variables para almacenar datos
energias = [];
gradientes = [];  % Norma del gradiente total
coordenadas = []; % Coordenadas atómicas [Ox,Oy,Oz,H1x,...]
pasos_optim = 0;

% Contadores para depuración
n_energias = 0;
n_coords = 0;
n_grads = 0;

fprintf("Analizando %d líneas del archivo...\n", length(lines));

% Bucle principal
for i = 1:length(lines)
    line = lines{i};

    % === 1. Extraer energía DFT ===
    % Patrón flexible para "Total DFT energy = -76.xxxx Ha"
    idx = strfind(line, "Total DFT energy =");
    if ~isempty(idx)
        % Extraer número después de "energy ="
        tokens = strsplit(line);
        for k = 1:length(tokens)
            if strcmp(tokens{k}, "energy") && k+2 <= length(tokens)
                energy_str = tokens{k+2};
                % Eliminar "Ha" si existe
                energy_str = strrep(energy_str, "Ha", "");
                energy_str = strrep(energy_str, ",", ".");  % Por si hay coma
                val = str2double(energy_str);
                if ~isnan(val)
                    energias(end+1) = val;
                    n_energias++;
                    fprintf(" Energía encontrada: %.8f Ha (línea %d)\n", val, i);
                endif
                break;
            endif
        endfor
    endif

    % === 2. Extraer coordenadas atómicas ===
    if ~isempty(strfind(line, "Output coordinates in angstroms"))
        fprintf("🔍 Coordenadas detectadas en línea %d\n", i);
        coords = zeros(3,3);
        found = 0;
        for j = i+2:i+10  % Las coordenadas vienen 2-3 líneas después
            if j > length(lines), break; endif
            coord_line = lines{j};
            parts = strsplit(strtrim(coord_line));
            if length(parts) >= 5
                % Intentar leer el número de átomo (1,2,3) y coordenadas
                atom_idx_str = parts{1};
                if isnumeric(str2double(atom_idx_str))
                    atom_idx = str2double(atom_idx_str);
                    if atom_idx >= 1 && atom_idx <= 3 && length(parts) >= 6
                        x = str2double(parts{4});
                        y = str2double(parts{5});
                        z = str2double(parts{6});
                        if ~isnan(x) && ~isnan(y) && ~isnan(z)
                            coords(atom_idx, :) = [x, y, z];
                            found++;
                            fprintf("   Átomo %d: (%.4f, %.4f, %.4f)\n", atom_idx, x, y, z);
                        endif
                    endif
                endif
            endif
        endfor
        if found == 3
            coordenadas(end+1, :) = coords(:)';  % Aplanar
            n_coords++;
            fprintf(" 3 coordenadas almacenadas.\n");
        else
            fprintf("  Solo %d coordenadas encontradas.\n", found);
        endif
    endif

    % === 3. Extraer gradientes ===
    if ~isempty(strfind(line, "DFT ENERGY GRADIENTS"))
        fprintf(" Gradientes detectados en línea %d\n", i);
        grads = zeros(3,3);
        found = 0;
        for j = i+3:i+10
            if j > length(lines), break; endif
            grad_line = lines{j};
            parts = strsplit(strtrim(grad_line));
            if length(parts) >= 7
                atom_idx_str = parts{1};
                if isnumeric(str2double(atom_idx_str))
                    atom_idx = str2double(atom_idx_str);
                    if atom_idx >= 1 && atom_idx <= 3
                        gx = str2double(parts{5});
                        gy = str2double(parts{6});
                        gz = str2double(parts{7});
                        if ~isnan(gx) && ~isnan(gy) && ~isnan(gz)
                            grads(atom_idx, :) = [gx, gy, gz];
                            found++;
                            fprintf("   Gradiente átomo %d: (%.6f, %.6f, %.6f)\n", atom_idx, gx, gy, gz);
                        endif
                    endif
                endif
            endif
        endfor
        if found == 3
            norm_grad = norm(grads(:));
            gradientes(end+1) = norm_grad;
            n_grads++;
            fprintf(" Gradientes almacenados. Norma total: %.3e Ha/Å\n", norm_grad);
        else
            fprintf("  Solo %d gradientes encontrados.\n", found);
        endif
    endif
endfor

% Mostrar resumen
fprintf("\n");
fprintf(" RESUMEN DE EXTRACCIÓN\n");
fprintf("Energías encontradas: %d\n", n_energias);
fprintf("Conjuntos de coordenadas: %d\n", size(coordenadas, 1));
fprintf("Conjuntos de gradientes: %d\n", length(gradientes));

if isempty(energias)
    warning("No se encontraron energías.");
else
    fprintf("Energía inicial: %.8f Ha\n", energias(1));
    fprintf("Energía final  : %.8f Ha\n", energias(end));
    fprintf("Disminución    : %.8f Ha\n", energias(1) - energias(end));
endif

% Calcular geometría si hay coordenadas
if ~isempty(coordenadas)
    n_steps = size(coordenadas, 1);
    distOH1 = zeros(n_steps, 1);
    distOH2 = zeros(n_steps, 1);
    angHOH = zeros(n_steps, 1);

    for i = 1:n_steps
        coord = reshape(coordenadas(i,:), 3, 3);
        O = coord(1, :);
        H1 = coord(2, :);
        H2 = coord(3, :);

        v1 = H1 - O;
        v2 = H2 - O;

        distOH1(i) = norm(v1);
        distOH2(i) = norm(v2);
        angHOH(i) = acosd(dot(v1, v2) / (norm(v1)*norm(v2)));
    endfor

    % Gráficas
    figure('Position', [100, 100, 800, 600]);

    subplot(2,2,1);
    plot(energias, 'o-');
    xlabel("Iteración");
    ylabel("Energía DFT (Ha)");
    title("Convergencia de energía");
    grid on;

    subplot(2,2,2);
    if ~isempty(gradientes)
        semilogy(gradientes, 's-');
        xlabel("Iteración");
        ylabel("||Gradient|| (Ha/Å)");
        title("Norma del gradiente");
        grid on;
    endif

    subplot(2,2,3);
    hold on;
    plot(distOH1, 'o-', 'DisplayName', 'O-H1');
    plot(distOH2, 's-', 'DisplayName', 'O-H2');
    xlabel("Paso");
    ylabel("Distancia (Å)");
    title("Longitudes de enlace");
    legend('show');
    grid on;
    hold off;

    subplot(2,2,4);
    plot(angHOH, '^-');
    xlabel("Paso");
    ylabel("Ángulo H-O-H (°)");
    title("Ángulo de enlace");
    grid on;

    % Resultados finales
    fprintf("\n--- Geometría final ---\n");
    fprintf("Distancia O-H1: %.4f Å\n", distOH1(end));
    fprintf("Distancia O-H2: %.4f Å\n", distOH2(end));
    fprintf("Ángulo H-O-H  : %.2f°\n", angHOH(end));
endif

% Guardar datos
data = struct();
data.energias = energias;
data.gradientes = gradientes;
data.coordenadas = coordenadas;
data.distOH1 = distOH1;
data.distOH2 = distOH2;
data.angHOH = angHOH;
save("datos_nwchem_extraidos.mat", "-struct", "data");

fprintf("\n Datos guardados en 'datos_nwchem_extraidos.mat'\n");
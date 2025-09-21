% parsear_nwchem.m
% Extracción de datos de salida de NWChem para H2O
% Compatible con Octave (todas las versiones)
% Detecta optimización DFT, energía, geometría, gradientes

clear; clc; close all;

% Nombre del archivo
filename = "frequences.nwo";

% Verificar archivo
if ~exist(filename, "file")
    error("No se encuentra el archivo: %s", filename);
endif

% Leer archivo línea por línea
fid = fopen(filename, "r");
if fid == -1
    error("No se pudo abrir el archivo.");
endif
lines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
lines = lines{1};  % Vector de cadenas

% Variables para almacenar datos
energias = [];
grad_norm = [];
distOH1 = [];
distOH2 = [];
angHOH = [];

% Bucle principal
for i = 1:length(lines)
    line = lines{i};

    % === 1. Extraer energía total DFT ===
    if ~isempty(strfind(line, "Total DFT energy ="))
        parts = strsplit(line);
        for k = 1:length(parts)
            if strcmp(parts{k}, "energy") && k+2 <= length(parts)
                energy_str = parts{k+2};
                energy_str = strrep(energy_str, "Ha", "");
                energy_str = strrep(energy_str, ",", ".");
                val = str2double(energy_str);
                if ~isnan(val)
                    energias(end+1) = val;
                    break;
                endif
            endif
        endfor
    endif

    % === 2. Extraer coordenadas atómicas ===
    if ~isempty(strfind(line, "Output coordinates in angstroms"))
        coords = zeros(3,3);
        found = 0;
        for j = i+2:i+8
            if j > length(lines), break; endif
            coord_line = lines{j};
            if ~isempty(strfind(coord_line, "Atomic Mass")), break; endif
            parts = strsplit(strtrim(coord_line));
            if length(parts) >= 6 && isnumeric(str2double(parts{1}))
                idx = str2double(parts{1});
                if idx >= 1 && idx <= 3
                    x = str2double(parts{4});
                    y = str2double(parts{5});
                    z = str2double(parts{6});
                    if ~isnan(x) && ~isnan(y) && ~isnan(z)
                        coords(idx,:) = [x,y,z];
                        found++;
                    endif
                endif
            endif
        endfor
        if found == 3
            O = coords(1,:); H1 = coords(2,:); H2 = coords(3,:);
            d1 = norm(H1 - O);
            d2 = norm(H2 - O);
            v1 = H1 - O; v2 = H2 - O;
            angle = acosd(dot(v1,v2)/(norm(v1)*norm(v2)));
            distOH1(end+1) = d1;
            distOH2(end+1) = d2;
            angHOH(end+1) = angle;
        endif
    endif

    % === 3. Extraer gradientes (DFT ENERGY GRADIENTS) ===
    if ~isempty(strfind(line, "DFT ENERGY GRADIENTS"))
        grads = zeros(3,3);
        found = 0;
        for j = i+3:i+10
            if j > length(lines), break; endif
            gline = lines{j};
            parts = strsplit(strtrim(gline));
            if length(parts) >= 7 && isnumeric(str2double(parts{1}))
                idx = str2double(parts{1});
                if idx >= 1 && idx <= 3
                    gx = str2double(parts{5});
                    gy = str2double(parts{6});
                    gz = str2double(parts{7});
                    if ~isnan(gx) && ~isnan(gy) && ~isnan(gz)
                        grads(idx,:) = [gx, gy, gz];
                        found++;
                    endif
                endif
            endif
        endfor
        if found == 3
            norm_g = norm(grads(:));
            grad_norm(end+1) = norm_g;
        endif
    endif
endfor

% Mostrar resultados
fprintf(" Extracción completada\n");
fprintf("Número de pasos de optimización: %d\n", length(energias));
if ~isempty(energias)
    fprintf("Energía inicial: %.8f Ha\n", energias(1));
    fprintf("Energía final  : %.8f Ha\n", energias(end));
    fprintf("ΔE             : %.2e Ha\n", energias(end) - energias(1));
endif

% Gráficas
figure('Position', [100, 100, 900, 600]);

subplot(2,3,1);
if ~isempty(energias)
    plot(energias, 'o-');
    xlabel("Paso");
    ylabel("Energía DFT (Ha)");
    title("Convergencia de energía");
    grid on;
endif

subplot(2,3,2);
if ~isempty(grad_norm)
    semilogy(grad_norm, 's-');
    xlabel("Paso");
    ylabel("||Gradient|| (Ha/Å)");
    title("Norma del gradiente");
    grid on;
endif

subplot(2,3,3);
if ~isempty(distOH1)
    hold on;
    plot(distOH1, 'o-', 'DisplayName', 'O-H1');
    plot(distOH2, 's-', 'DisplayName', 'O-H2');
    xlabel("Paso");
    ylabel("Distancia (Å)");
    title("Longitudes de enlace");
    legend('show');
    grid on;
    hold off;
endif

subplot(2,3,4);
if ~isempty(angHOH)
    plot(angHOH, '^-');
    xlabel("Paso");
    ylabel("Ángulo H-O-H (°)");
    title("Ángulo de enlace");
    ylim([80, 120]);
    grid on;
endif

% Tabla final
if ~isempty(distOH1) && ~isempty(angHOH)
    fprintf("\n--- Geometría final optimizada ---\n");
    fprintf("O-H1 = %.4f Å\n", distOH1(end));
    fprintf("O-H2 = %.4f Å\n", distOH2(end));
    fprintf("H-O-H = %.2f°\n", angHOH(end));
endif

% Guardar datos
data = struct();
data.energias = energias;
data.grad_norm = grad_norm;
data.distOH1 = distOH1;
data.distOH2 = distOH2;
data.angHOH = angHOH;
save("datos_nwchem_extraidos.mat", "-struct", "data");

fprintf("\n Datos guardados en 'datos_nwchem_extraidos.mat'\n");

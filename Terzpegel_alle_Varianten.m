%% ============================================================
%  Terzpegel-Berechnung aus RIR-Dateien
%  Ausgabe: Terzpegel-Matrix [14 x 48]
%  Export: Excel (.xlsx)
% ============================================================

clear;
clc;

%% ---------------- Einstellungen ----------------
nMess = 14;
fs = 500e3;        % Abtastrate [Hz]

% Normgerechte 1/3-Oktav-Mittenfrequenzen (IEC 61260)
f_terz = double([ ...
    560 625 710 800 900 1000 1120 1250 1400 1600 1800 2000 ...
    2240 2500 2800 3150 3550 4000 4500 5000 5600 6300 7100 8000 ...
    9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 ...
    28000 31500 35500 40000 45000 50000 56000 63000 71000 80000 ...
    90000 100000 112000 126000 ]);

filterOrder = 4;
ref = 1;

% Alle Variante_x Ordner finden
variantFolders = dir('Variante_*');
variantFolders = variantFolders([variantFolders.isdir]);

% Nur Ordner ab 'Variante_1' verarbeiten
startFolder = 'Variante_1_neu';
variantNames = {variantFolders.name};
startIdx = find(strcmp(variantNames, startFolder), 1);

if isempty(startIdx)
    warning('Startordner "%s" nicht gefunden. Alle Varianten werden verarbeitet.', startFolder);
else
    variantFolders = variantFolders(startIdx:end);
    fprintf('Starte Verarbeitung ab: %s\n', startFolder);
end

if isempty(variantFolders)
    error('Keine Variante_x Ordner gefunden!');
end

fprintf('Gefundene Varianten: %d\n', length(variantFolders));
for v = 1:length(variantFolders)
    fprintf('  - %s\n', variantFolders(v).name);
end

%% ---------------- Verarbeitung aller Varianten ----------------

for varIdx = 1:length(variantFolders)

    variantName = variantFolders(varIdx).name;
    fprintf('\n========================================\n');
    fprintf('Verarbeite %s...\n', variantName);
    fprintf('========================================\n');

    % Dateinamenmuster für diese Variante
    filePattern = sprintf('%s/%s,Pos_%%d.mat', variantName, variantName);

%% ---------------- Speicher ----------------
nFreq = numel(f_terz);
L = NaN(nMess, nFreq);   % [Messung x Terzband]

%% ---------------- Verarbeitung ----------------
for i = 1:nMess

    filename = sprintf(filePattern, i);

    if ~exist(filename, 'file')
        warning('Datei nicht gefunden: %s', filename);
        continue;
    end

    data = load(filename, 'RIR');

    if ~isfield(data, 'RIR')
        error('Datei %s enthält keinen Vektor "RIR".', filename);
    end

    ir = data.RIR(:);   % Spaltenvektor

    for k = 1:nFreq
        % --- Terzband-Grenzen ---
        f1 = f_terz(k) / 2^(1/6);
        f2 = f_terz(k) * 2^(1/6);

        Wn = [f1 f2] / (fs/2);

        if Wn(2) >= 1
            continue;   % oberhalb Nyquist
        end

        % --- Terzbandfilter ---
        [b,a] = butter(filterOrder/2, Wn, 'bandpass');
        ir_filt = filtfilt(b, a, ir);

        % --- Terzpegel berechnen ---
        L(i,k) = 10*log10(sum(ir_filt.^2)/ref^2 + eps);
    end
end

%% ---------------- Excel-Export ----------------

% Zeilen- und Spaltennamen (kurz & MATLAB-sicher)
rowNames = arrayfun(@(x) sprintf('Pos_%02d', x), 1:nMess, ...
                    'UniformOutput', false);

colNames = arrayfun(@(f) sprintf('F%.0f', f), f_terz, ...
                    'UniformOutput', false);

T_Terzpegel = array2table(L, ...
    'RowNames', rowNames, ...
    'VariableNames', colNames);

% Excel in Varianten-Ordner speichern
excelFile = fullfile(variantName, 'Terzpegel_dBFS.xlsx');
writetable(T_Terzpegel, excelFile, 'WriteRowNames', true);

fprintf('Excel-Datei erstellt: %s\n', excelFile);

%% ---------------- Plot: Terzpegel pro Position ----------------

% Globale Y-Achsen-Grenzen bestimmen (für Vergleichbarkeit)
y_min = min(L(:), [], 'omitnan');
y_max = max(L(:), [], 'omitnan');
y_range = [floor(y_min/10)*10, ceil(y_max/10)*10];

for i = 1:nMess
    fig = figure('Position', [100, 100, 1000, 500]);

    % Terzpegel als Stairs-Diagramm
    stairs(f_terz, L(i,:), ...
           'LineWidth', 2, ...
           'Color', [0 0.4470 0.7410]);

    grid on;
    set(gca, 'XScale', 'log');
    xlabel('Frequenz [Hz]');
    ylabel('Terzpegel [dBFS]');
    title(sprintf('Terzpegel - Position %02d', i));
    xlim([min(f_terz) max(f_terz)]);
    ylim(y_range);

    % Formatierung der x-Achse
    xticks([500 1000 2000 5000 10000 20000 50000 100000]);
    xticklabels({'500', '1k', '2k', '5k', '10k', '20k', '50k', '100k'});

    % Plot in Varianten-Ordner speichern
    filename = fullfile(variantName, sprintf('Terzpegel_Pos_%02d', i));
    saveas(fig, [filename '.png']);
    saveas(fig, [filename '.fig']);

    fprintf('Plot erstellt: %s.png\n', filename);
end

% Alle Figuren schließen
close all;

disp('Alle Terzpegel-Plots erstellt und gespeichert.');

%% ---------------- Mittelwert ausgewählter Positionen ----------------

% Positionen die gemittelt werden sollen
positions_to_average = [5, 6, 7, 9, 10, 11, 13, 14];

% Mittelwert über die ausgewählten Positionen berechnen
L_mean = mean(L(positions_to_average, :), 1, 'omitnan');

% Plot erstellen
fig_mean = figure('Position', [100, 100, 1000, 500]);

stairs(f_terz, L_mean, ...
       'LineWidth', 2.5, ...
       'Color', [0.8500 0.3250 0.0980]);  % Orange für Mittelwert

grid on;
set(gca, 'XScale', 'log');
xlabel('Frequenz [Hz]');
ylabel('Terzpegel [dBFS]');
title(sprintf('%s - Gemittelter Terzpegel (Pos. 5,6,7,9,10,11,13,14)', variantName));
xlim([min(f_terz) max(f_terz)]);
ylim(y_range);

% Formatierung der x-Achse
xticks([500 1000 2000 5000 10000 20000 50000 100000]);
xticklabels({'500', '1k', '2k', '5k', '10k', '20k', '50k', '100k'});

% Plot in Varianten-Ordner speichern
filename_mean = fullfile(variantName, 'Terzpegel_Mittelwert');
saveas(fig_mean, [filename_mean '.png']);
saveas(fig_mean, [filename_mean '.fig']);

fprintf('Mittelwert-Plot erstellt: %s.png\n', filename_mean);
close(fig_mean);

fprintf('\n%s abgeschlossen!\n', variantName);

end  % Ende der Varianten-Schleife

disp('========================================');
disp('Alle Varianten erfolgreich verarbeitet!');
disp('========================================');

%% ============================================================
%  Terzband-Nachhallzeit (RT60) aus RIR-Dateien
%  Methode: Schroeder-Integration + T20
%  Ausgabe: RT60-Matrix [14 x 48]
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

% Alle Variante_x Ordner finden
variantFolders = dir('Variante_*');
variantFolders = variantFolders([variantFolders.isdir]);

% Nur Ordner ab 'Variante_X' verarbeiten
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
RT60_data = NaN(nMess, nFreq);   % [Messung x Terzband]

%% ---------------- Verarbeitung ----------------
for i = 1:nMess

    filename = sprintf(filePattern, i);
    data = load(filename,'RIR');

    if ~isfield(data,'RIR')
        error('Datei %s enthält keinen Vektor "RIR".', filename);
    end

    ir = data.RIR(:);

    % --- Lundeby-Truncation: Finde Rauschgrenze ---
    % Einfache Methode: Finde wo Signal unter Rauschpegel fällt
    ir_abs = abs(ir);
    max_amp = max(ir_abs);

    % Rauschpegel aus letzten 10% schätzen
    noise_samples = ir_abs(end-round(length(ir)*0.1):end);
    noise_level = mean(noise_samples) + 3*std(noise_samples);

    % Finde letzten signifikanten Peak (mindestens 10x über Rauschpegel)
    threshold = max(noise_level * 10, max_amp * 0.001);  % mindestens -60 dB
    sig_idx = find(ir_abs > threshold, 1, 'last');

    if isempty(sig_idx) || sig_idx < 100
        sig_idx = length(ir);  % Falls keine Truncation möglich
    else
        % Sicherheitsmarge: 20% zusätzlich
        sig_idx = min(length(ir), round(sig_idx * 1.2));
    end

    % Truncierte Impulsantwort
    ir_trunc = ir(1:sig_idx);
    t = (0:length(ir_trunc)-1)'/fs;

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
        ir_filt = filtfilt(b,a,ir_trunc);

        % --- Schroeder-Integration ---
        E = flipud(cumsum(flipud(ir_filt.^2)));

        % Energie muss > 0 sein für Log
        if max(E) == 0
            continue;
        end

        E = E / max(E);  % Normierung
        E_dB = 10*log10(E + eps);

        % --- T20-Bereich ---
        idx = find(E_dB <= -5 & E_dB >= -25);

        if numel(idx) < 10
            continue;
        end

        % --- Lineare Regression ---
        p = polyfit(t(idx), E_dB(idx), 1);

        % Steigung muss negativ sein (Abfall)
        if p(1) >= 0
            continue;
        end

        % --- RT60 ---
        RT60_data(i,k) = -60 / p(1);   % Sekunden
    end
end

%% ---------------- Excel-Export ----------------

% Zeilen- und Spaltennamen (kurz & MATLAB-sicher)
rowNames = arrayfun(@(x) sprintf('Pos_%02d', x), 1:nMess, ...
                    'UniformOutput', false);

colNames = arrayfun(@(f) sprintf('F%.0f', f), f_terz, ...
                    'UniformOutput', false);

T_RT60 = array2table(RT60_data, ...
    'RowNames', rowNames, ...
    'VariableNames', colNames);

% Excel in Varianten-Ordner speichern
excelFile = fullfile(variantName, 'RT60_Terzband.xlsx');
writetable(T_RT60, excelFile, 'WriteRowNames', true);

%% ---------------- Plot: RT60 pro Position ----------------

% Globale Y-Achsen-Grenzen bestimmen (für Vergleichbarkeit)
y_min_RT60 = min(RT60_data(:), [], 'omitnan');
y_max_RT60 = max(RT60_data(:), [], 'omitnan');
y_range_RT60 = [floor(y_min_RT60*10)/10, ceil(y_max_RT60*10)/10];

for i = 1:nMess
    fig = figure('Position', [100, 100, 1000, 500]);

    % RT60 als Stairs-Diagramm
    stairs(f_terz, RT60_data(i,:), ...
           'LineWidth', 2, ...
           'Color', [0 0.4470 0.7410]);

    grid on;
    set(gca, 'XScale', 'log');
    xlabel('Frequenz [Hz]');
    ylabel('Nachhallzeit RT60 [s]');
    title(sprintf('RT60 - Position %02d', i));
    xlim([min(f_terz) max(f_terz)]);
    ylim(y_range_RT60);

    % Formatierung der x-Achse
    xticks([500 1000 2000 5000 10000 20000 50000 100000]);
    xticklabels({'500', '1k', '2k', '5k', '10k', '20k', '50k', '100k'});

    % Plot in Varianten-Ordner speichern
    filename = fullfile(variantName, sprintf('RT60_Pos_%02d', i));
    saveas(fig, [filename '.png']);
    saveas(fig, [filename '.fig']);

    fprintf('Plot erstellt: %s.png\n', filename);
end

% Alle RT60-Figuren schließen
close all;

disp('Alle RT60-Plots erstellt und gespeichert.');

%% ---------------- Mittelwert ausgewählter Positionen ----------------

% Positionen die gemittelt werden sollen
positions_to_average = [5, 6, 7, 9, 10, 11, 13, 14];

% Mittelwert über die ausgewählten Positionen berechnen
RT60_mean = mean(RT60_data(positions_to_average, :), 1, 'omitnan');

% Plot erstellen
fig_mean = figure('Position', [100, 100, 1000, 500]);

stairs(f_terz, RT60_mean, ...
       'LineWidth', 2.5, ...
       'Color', [0.8500 0.3250 0.0980]);  % Orange für Mittelwert

grid on;
set(gca, 'XScale', 'log');
xlabel('Frequenz [Hz]');
ylabel('Nachhallzeit RT60 [s]');
title(sprintf('%s - Gemittelte RT60 (Pos. 5,6,7,9,10,11,13,14)', variantName));
xlim([min(f_terz) max(f_terz)]);
ylim(y_range_RT60);

% Formatierung der x-Achse
xticks([500 1000 2000 5000 10000 20000 50000 100000]);
xticklabels({'500', '1k', '2k', '5k', '10k', '20k', '50k', '100k'});

% Plot in Varianten-Ordner speichern
filename_mean = fullfile(variantName, 'RT60_Mittelwert');
saveas(fig_mean, [filename_mean '.png']);
saveas(fig_mean, [filename_mean '.fig']);

fprintf('Mittelwert-Plot erstellt: %s.png\n', filename_mean);
close(fig_mean);

fprintf('\n%s abgeschlossen!\n', variantName);

end  % Ende der Varianten-Schleife

disp('========================================');
disp('Alle Varianten erfolgreich verarbeitet!');
disp('========================================');

%% ---------------- Vergleich der Mittelwerte aller Varianten ----------------

% Varianten die verglichen werden sollen
variantsToCompare = {'Variante_1_neu', 'Variante_2', 'Variante_3', 'Variante_4'};
nVariants = length(variantsToCompare);

% Positionen die gemittelt werden sollen (wie oben definiert)
positions_to_average = [5, 6, 7, 9, 10, 11, 13, 14];

% Speicher für Mittelwerte aller Varianten
RT60_means = NaN(nVariants, nFreq);

% Lade die gespeicherten RT60-Daten aus den Excel-Dateien
fprintf('\n========================================\n');
fprintf('Lade RT60-Daten für Variantenvergleich...\n');
fprintf('========================================\n');

for v = 1:nVariants
    varName = variantsToCompare{v};
    excelPath = fullfile(varName, 'RT60_Terzband.xlsx');

    if exist(excelPath, 'file')
        % Excel-Datei einlesen (erste Spalte sind Zeilennamen)
        T = readtable(excelPath, 'ReadRowNames', true);
        RT60_matrix = table2array(T);

        % Mittelwert über ausgewählte Positionen berechnen
        RT60_means(v, :) = mean(RT60_matrix(positions_to_average, :), 1, 'omitnan');

        fprintf('✓ %s geladen\n', varName);
    else
        warning('Datei nicht gefunden: %s', excelPath);
    end
end

% Plot erstellen: Vergleich aller Varianten
fig_compare = figure('Position', [100, 100, 1200, 600]);

% Farben für verschiedene Varianten
colors = [
    0.0000 0.4470 0.7410;  % Blau - Variante 1_neu
    0.8500 0.3250 0.0980;  % Orange - Variante 2
    0.9290 0.6940 0.1250;  % Gelb - Variante 3
    0.4940 0.1840 0.5560   % Lila - Variante 4
];

hold on;
legend_entries = {};

for v = 1:nVariants
    stairs(f_terz, RT60_means(v, :), ...
           'LineWidth', 2.5, ...
           'Color', colors(v, :), ...
           'DisplayName', variantsToCompare{v});
    legend_entries{v} = strrep(variantsToCompare{v}, '_', ' ');
end

hold off;

grid on;
set(gca, 'XScale', 'log');
xlabel('Frequenz [Hz]', 'FontSize', 12);
ylabel('Nachhallzeit RT60 [s]', 'FontSize', 12);
title('Vergleich RT60-Mittelwerte aller Varianten (Pos. 5,6,7,9,10,11,13,14)', 'FontSize', 14);
xlim([min(f_terz) max(f_terz)]);

% Y-Achse: Dynamisch an Daten anpassen
y_min_all = min(RT60_means(:), [], 'omitnan');
y_max_all = max(RT60_means(:), [], 'omitnan');
ylim([floor(y_min_all*10)/10, ceil(y_max_all*10)/10]);

% Formatierung der x-Achse
xticks([500 1000 2000 5000 10000 20000 50000 100000]);
xticklabels({'500', '1k', '2k', '5k', '10k', '20k', '50k', '100k'});

% Legende
legend(legend_entries, 'Location', 'best', 'FontSize', 11);

% Plot speichern
filename_compare = 'Vergleich_RT60_Mittelwerte_alle_Varianten';
saveas(fig_compare, [filename_compare '.png']);
saveas(fig_compare, [filename_compare '.fig']);

fprintf('\n========================================\n');
fprintf('Vergleichsplot erstellt: %s.png\n', filename_compare);
fprintf('========================================\n');

close(fig_compare);


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---

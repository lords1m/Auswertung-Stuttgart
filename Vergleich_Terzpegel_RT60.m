%% ============================================================
%  Vergleich: Terzpegel & RT60 für Position 7
%  Varianten: 1, 2, 4
%  Ausgabe: Kombinierte Plots mit zwei Y-Achsen
% ============================================================

clear;
clc;

%% ---------------- Einstellungen ----------------

% Zu vergleichende Varianten
varianten = {'Variante_1_neu', 'Variante_3', 'Variante_4'};
position = 9;  % Position die verglichen werden soll


f_terz = double([ ...
    560 625 710 800 900 1000 1120 1250 1400 1600 1800 2000 ...
    2240 2500 2800 3150 3550 4000 4500 5000 5600 6300 7100 8000 ...
    9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 ...
    28000 31500 35500 40000 45000 50000 56000 63000 71000 80000 ...
    90000 100000 112000 126000 ]);

% Farben für die Varianten
colors = [0 0.4470 0.7410;      % Blau
          0.8500 0.3250 0.0980;  % Orange
          0.9290 0.6940 0.1250]; % Gelb

%% ---------------- Daten laden ----------------

nVarianten = length(varianten);
nFreq = length(f_terz);

% Speicher für Daten
Terzpegel_data = NaN(nVarianten, nFreq);
RT60_data = NaN(nVarianten, nFreq);

for v = 1:nVarianten
    variantName = varianten{v};

    % Terzpegel laden
    terzpegel_file = fullfile(variantName, 'Terzpegel_dBFS.xlsx');
    if exist(terzpegel_file, 'file')
        T_terz = readtable(terzpegel_file, 'ReadRowNames', true);
        row_name = sprintf('Pos_%02d', position);
        if ismember(row_name, T_terz.Properties.RowNames)
            Terzpegel_data(v, :) = table2array(T_terz(row_name, :));
            fprintf('✓ Terzpegel geladen: %s, Position %d\n', variantName, position);
        else
            warning('Position %d nicht gefunden in %s', position, terzpegel_file);
        end
    else
        warning('Datei nicht gefunden: %s', terzpegel_file);
    end

    % RT60 laden
    rt60_file = fullfile(variantName, 'RT60_Terzband.xlsx');
    if exist(rt60_file, 'file')
        T_rt60 = readtable(rt60_file, 'ReadRowNames', true);
        row_name = sprintf('Pos_%02d', position);
        if ismember(row_name, T_rt60.Properties.RowNames)
            RT60_data(v, :) = table2array(T_rt60(row_name, :));
            fprintf('✓ RT60 geladen: %s, Position %d\n', variantName, position);
        else
            warning('Position %d nicht gefunden in %s', position, rt60_file);
        end
    else
        warning('Datei nicht gefunden: %s', rt60_file);
    end
end

%% ---------------- Plot 1: Separater Plot für Terzpegel und RT60 ----------------

fig1 = figure('Position', [100, 100, 1200, 800]);

% Subplot 1: Terzpegel
subplot(2, 1, 1);
hold on;
grid on;

for v = 1:nVarianten
    stairs(f_terz, Terzpegel_data(v, :), ...
           'LineWidth', 2, ...
           'Color', colors(v, :), ...
           'DisplayName', varianten{v});
end

set(gca, 'XScale', 'log');
xlabel('Frequenz [Hz]');
ylabel('Terzpegel [dBFS]');
title(sprintf('Terzpegel-Vergleich - Position %d', position));
xlim([min(f_terz) max(f_terz)]);
xticks([500 1000 2000 5000 10000 20000 50000 100000]);
xticklabels({'500', '1k', '2k', '5k', '10k', '20k', '50k', '100k'});
legend('Location', 'best');
hold off;

% Subplot 2: RT60
subplot(2, 1, 2);
hold on;
grid on;

for v = 1:nVarianten
    stairs(f_terz, RT60_data(v, :), ...
           'LineWidth', 2, ...
           'Color', colors(v, :), ...
           'DisplayName', varianten{v});
end

set(gca, 'XScale', 'log');
xlabel('Frequenz [Hz]');
ylabel('Nachhallzeit RT60 [s]');
title(sprintf('RT60-Vergleich - Position %d', position));
xlim([min(f_terz) max(f_terz)]);
xticks([500 1000 2000 5000 10000 20000 50000 100000]);
xticklabels({'500', '1k', '2k', '5k', '10k', '20k', '50k', '100k'});
legend('Location', 'best');
hold off;

% Speichern - Dateiname erstellen mit Varianten und Position
varianten_str = strjoin(strrep(varianten, 'Variante_', 'V'), '_');
filename_base = sprintf('Vergleich_Terzpegel_RT60_%s_Pos%02d', varianten_str, position);
saveas(fig1, [filename_base '.png']);
saveas(fig1, [filename_base '.fig']);
fprintf('\nPlot 1 gespeichert: %s.png\n', filename_base);



%% ---------------- Zusammenfassung ----------------

fprintf('\n========================================\n');
fprintf('Vergleichsanalyse abgeschlossen!\n');
fprintf('========================================\n');
fprintf('Position: %d\n', position);
fprintf('Varianten: %s\n', strjoin(varianten, ', '));
fprintf('\nErstellt:\n');
fprintf('  - %s.png (2 Subplots)\n', filename_base);
fprintf('========================================\n');

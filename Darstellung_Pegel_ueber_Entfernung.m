%% ============================================================
%  Darstellung: Schallpegel über Entfernung von der Quelle
%  Ergebnis: Plot L vs. s mit idealer 1/r-Linie zum Vergleich
%
%  Zeigt den Summenpegel (oder Terzpegel bei gewählter Frequenz)
%  als Funktion der Entfernung mit theoretischer Freifeld-Dämpfung
% ============================================================

clear;
clc;
close all;

% Stelle sicher, dass wir im richtigen Verzeichnis sind
scriptDir = fileparts(mfilename('fullpath'));
if ~isempty(scriptDir)
    cd(scriptDir);
end
fprintf('Arbeitsverzeichnis: %s\n', pwd);

%% ---------------- Einstellungen (anpassen) ----------------

% --- Varianten auswählen ---
% Kann eine oder mehrere Varianten sein
varianten = {'Variante_1_neu', 'Variante_2'};
% varianten = {'Variante_1_neu'}; % Einzelne Variante

% --- Messpositionen und Entfernungen ---
% Die Quellenposition ist Q1 (Position 1 in Reihe 4)
% Entfernungen in Metern (bitte anpassen basierend auf tatsächlichem Setup)

% Beispiel-Layout (4x4 Grid):
% Annahme: Quadratisches Gitter mit 1m Abstand zwischen Positionen
% Q1 ist an Position (0, 0) in Reihe 4, Spalte 1

positions_info = struct();

% Reihe 1 (y = 3m)
positions_info(1).pos = 1;  positions_info(1).x = 0;   positions_info(1).y = 1.2;
positions_info(2).pos = 2;  positions_info(2).x = 0.3;   positions_info(2).y = 1.2;
positions_info(3).pos = 3;  positions_info(3).x = 0.6;   positions_info(3).y = 1.2;
positions_info(4).pos = 4;  positions_info(4).x = 1.2;   positions_info(4).y = 1.2;

% Reihe 2 (y = 2m)
positions_info(5).pos = 5;  positions_info(5).x = 0;   positions_info(5).y = 0.6;
positions_info(6).pos = 6;  positions_info(6).x = 0.3;   positions_info(6).y = 0.6;
positions_info(7).pos = 7;  positions_info(7).x = 0.6;   positions_info(7).y = 0.6;
positions_info(8).pos = 8;  positions_info(8).x = 1.2;   positions_info(8).y = 0.6;

% Reihe 3 (y = 1m)
positions_info(9).pos = 9;   positions_info(9).x = 0;   positions_info(9).y = 0.3;
positions_info(10).pos = 10; positions_info(10).x = 0.3;  positions_info(10).y = 0.3;
positions_info(11).pos = 11; positions_info(11).x = 0.6;  positions_info(11).y = 0.3;
positions_info(12).pos = 12; positions_info(12).x = 1.2;  positions_info(12).y = 0.3;

% Reihe 4 (y = 0m) - Q1 ist die Quelle bei (0,0)
positions_info(13).pos = 13; positions_info(13).x = 0.3;  positions_info(13).y = 0;
positions_info(14).pos = 14; positions_info(14).x = 0.6;  positions_info(14).y = 0;
positions_info(15).pos = 15; positions_info(15).x = 1.2;  positions_info(15).y = 0;

% Quellenposition (Position Q1 = Pos 1 in den Dateien mit "Q")
source_x = 0;
source_y = 0;

% Berechne Entfernungen von der Quelle
for i = 1:length(positions_info)
    dx = positions_info(i).x - source_x;
    dy = positions_info(i).y - source_y;
    positions_info(i).distance = sqrt(dx^2 + dy^2);
end

% --- Analyse-Modus ---
% 'summenpegel': Gesamtpegel über alle Frequenzen
% 'terzband': Pegel bei einer spezifischen Terzbandfrequenz
analyse_modus = 'summenpegel';
% analyse_modus = 'terzband';

% Wenn Terzband-Modus gewählt, welche Frequenz?
target_frequency = 16000; % Hz (z.B. 16 kHz)

% --- Darstellungs-Modus ---
% 'pegel': Darstellung in dBFS (logarithmisch)
% 'energie': Darstellung als relative Intensität I/I_ref (linear)
darstellung_modus = 'energie';
% darstellung_modus = 'energie';

% --- Referenzpegel für 1/r-Kurve ---
% Der Pegel an der Referenzentfernung (wird automatisch berechnet)
% Alternativ manuell setzen für Vergleich
use_auto_reference = true; % Automatisch: Nimmt gemessenen Pegel bei kleinster Entfernung
manual_reference_level = 0; % dBFS bei Referenzentfernung (nur wenn use_auto_reference = false)
manual_reference_distance = 0.5; % m (nur wenn use_auto_reference = false)

% --- Allgemeine Einstellungen ---
dataDir = 'data'; % Datenordner mit allen .mat-Dateien
fs = 500e3; % 500 kHz Abtastrate

% Normgerechte 1/3-Oktav-Mittenfrequenzen (IEC 61260)
% Begrenzt auf 4 kHz - 60 kHz
f_terz = double([ ...
    4000 4500 5000 5600 6300 7100 8000 ...
    9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 ...
    28000 31500 35500 40000 45000 50000 56000 ]);

%% ---------------- Prüfe data-Ordner ----------------
if ~exist(dataDir, 'dir')
    error('Datenordner "%s" nicht gefunden!', dataDir);
end

%% ---------------- Dateien einlesen und globale Referenz finden ----------------
dirInfo = dir(fullfile(dataDir, '*.mat'));
matFiles = {dirInfo.name};
if isempty(matFiles)
    error('Keine .mat-Dateien im Ordner "%s" gefunden!', dataDir);
end

FS_global = 0;
fprintf('Suche globale Referenz (FS) in allen .mat Dateien...\n');
for i = 1:numel(matFiles)
    try
        S = load(fullfile(dataDir, matFiles{i}));
        ir = extractIR(S);
        if ~isempty(ir)
            [ir_trunc, ~, ~, ~, ~, ~] = truncateIR(ir);
            N = numel(ir_trunc);
            H_mag = abs(fft(ir_trunc, N));
            FS_global = max(FS_global, max(H_mag(1:floor(N/2)+1)));
        end
    catch ME
        warning('Datei %s konnte nicht für FS-Bestimmung gelesen werden: %s', matFiles{i}, ME.message);
    end
end

if FS_global == 0
    error('Keine gültigen Impulsantworten gefunden. FS_Global konnte nicht bestimmt werden.');
end
fprintf('✓ Globale Referenz (FS) ermittelt: %g\n', FS_global);

%% ---------------- Pegel für alle Positionen berechnen ----------------
fprintf('\n=== Berechne Pegel für alle Varianten und Positionen ===\n');

% Struktur zum Speichern der Ergebnisse
results = struct();

for v_idx = 1:numel(varianten)
    variante = varianten{v_idx};
    fprintf('\n--- Variante: %s ---\n', variante);

    results(v_idx).variante = variante;
    results(v_idx).positions = [];
    results(v_idx).distances = [];
    results(v_idx).levels = [];
    results(v_idx).energies = []; % Energie aus Impulsantwort

    for p_idx = 1:length(positions_info)
        position = positions_info(p_idx).pos;
        distance = positions_info(p_idx).distance;

        if distance == 0
            % Quellenposition überspringen
            continue;
        end

        fprintf('  Position %d (Entfernung: %.2f m)... ', position, distance);

        % Berechne Pegel
        if strcmp(analyse_modus, 'summenpegel')
            level = getSummenpegel(variante, position, dataDir, matFiles, fs, f_terz, FS_global);
        else
            % Terzband-Modus
            terzpegel = getTerzpegel(variante, position, dataDir, matFiles, fs, f_terz, FS_global);
            [~, freq_idx] = min(abs(f_terz - target_frequency));
            level = terzpegel(freq_idx);
        end

        if ~isnan(level)
            % Hole Impulsantwort für Energieberechnung (wird für Energie-Modus benötigt)
            energy = NaN;

            % Suche nach der Datei mit verschiedenen Namensmustern (wie getTerzpegel)
            posStr = num2str(position);
            patterns = {
                sprintf('%s,Pos_%s.mat', variante, posStr),
                sprintf('%s_neu,Pos_%s.mat', variante, posStr),
                sprintf('%s,Pos%s.mat', variante, posStr)
            };

            filepath = '';
            for p = 1:length(patterns)
                for f = 1:length(matFiles)
                    if strcmpi(matFiles{f}, patterns{p})
                        filepath = fullfile(dataDir, matFiles{f});
                        break;
                    end
                end
                if ~isempty(filepath)
                    break;
                end
            end

            if ~isempty(filepath)
                try
                    S = load(filepath);
                    ir = extractIR(S);
                    if ~isempty(ir)
                        [ir_trunc, ~, ~, ~, ~, ~] = truncateIR(ir);
                        % Berechne Energie als Summe der quadrierten Samples
                        energy = sum(ir_trunc.^2);
                    end
                catch ME
                    % Fehler beim Laden, energy bleibt NaN
                    % fprintf('Fehler beim Laden von %s: %s\n', filepath, ME.message);
                end
            end

            results(v_idx).positions(end+1) = position;
            results(v_idx).distances(end+1) = distance;
            results(v_idx).levels(end+1) = level;
            results(v_idx).energies(end+1) = energy;

            if ~isnan(energy)
                fprintf('✓ L = %.2f dBFS, E = %.2e\n', level, energy);
            else
                fprintf('✓ L = %.2f dBFS (keine Energie-Daten)\n', level);
            end
        else
            fprintf('✗ Keine Daten\n');
        end
    end

    if isempty(results(v_idx).distances)
        warning('Keine gültigen Daten für Variante %s gefunden!', variante);
    end
end

%% ---------------- Ideale 1/r-Kurve berechnen ----------------
fprintf('\n=== Berechne ideale 1/r-Dämpfungskurve ===\n');

% Finde minimale und maximale Entfernung über alle Messungen
all_distances = [];
all_levels = [];
for v_idx = 1:numel(results)
    all_distances = [all_distances, results(v_idx).distances];
    all_levels = [all_levels, results(v_idx).levels];
end

if isempty(all_distances)
    error('Keine gültigen Messungen gefunden. Abbruch.');
end

min_dist = min(all_distances);
max_dist = max(all_distances);

% Erstelle Entfernungsvektor für ideale Kurve
distance_ideal = linspace(min_dist, max_dist, 200);

% Berechne Referenzpegel
if use_auto_reference
    % Nimm gemessenen Pegel bei kleinster Entfernung
    % (Mittelwert über alle Varianten bei dieser Entfernung)
    levels_at_min_dist = [];
    for v_idx = 1:numel(results)
        idx = find(abs(results(v_idx).distances - min_dist) < 0.01);
        if ~isempty(idx)
            levels_at_min_dist = [levels_at_min_dist, results(v_idx).levels(idx(1))];
        end
    end
    if ~isempty(levels_at_min_dist)
        L_ref = mean(levels_at_min_dist);
        r_ref = min_dist;
        fprintf('Automatische Referenz: L_ref = %.2f dBFS bei r_ref = %.2f m\n', L_ref, r_ref);
    else
        error('Keine Messung bei minimaler Entfernung gefunden.');
    end
else
    L_ref = manual_reference_level;
    r_ref = manual_reference_distance;
    fprintf('Manuelle Referenz: L_ref = %.2f dBFS bei r_ref = %.2f m\n', L_ref, r_ref);
end

% Ideale 1/r-Dämpfung für Halbraum: L(r) = L_ref - 20*log10(r/r_ref)
% Im Halbraum nimmt der Schalldruckpegel um 6 dB pro Abstandsverdopplung ab
L_ideal = L_ref - 20*log10(distance_ideal / r_ref);

%% ---------------- Umrechnung in Energie-Darstellung (falls gewählt) ----------------
if strcmp(darstellung_modus, 'energie')
    fprintf('\n=== Wechsel zu Energie-Darstellung ===\n');

    % Verwende die berechneten Energien aus den Impulsantworten
    % Berechne Referenzenergie (Energie bei r_ref)
    E_ref = NaN;
    for v_idx = 1:numel(results)
        idx = find(abs(results(v_idx).distances - r_ref) < 0.01);
        if ~isempty(idx)
            E_ref = results(v_idx).energies(idx(1));
            break;
        end
    end

    % Falls keine Energie exakt bei r_ref, nimm die bei minimaler Entfernung
    if isnan(E_ref)
        for v_idx = 1:numel(results)
            idx = find(abs(results(v_idx).distances - min_dist) < 0.01);
            if ~isempty(idx)
                E_ref = results(v_idx).energies(idx(1));
                break;
            end
        end
    end

    if isnan(E_ref)
        warning('Keine Referenzenergie gefunden! Alle energy-Werte sind NaN.');
        warning('Wechsle zurück zu Pegel-Darstellung.');
        darstellung_modus = 'pegel';
        % Neuberechnung nicht nötig, da levels bereits Pegel enthalten
    end
end

if strcmp(darstellung_modus, 'energie')
    % Nur wenn Energie-Modus erfolgreich war

    fprintf('  Referenzenergie E_ref = %.2e bei r_ref = %.2f m\n', E_ref, r_ref);

    % Ersetze levels durch relative Energien E/E_ref
    for v_idx = 1:numel(results)
        if ~isempty(results(v_idx).energies)
            results(v_idx).levels = results(v_idx).energies / E_ref;
        end
    end

    % Berechne ideale Energiekurve für Halbraum: E(r) = E_ref * (r_ref/r)^2
    % (Im Halbraum: Intensität ~ 1/r², Energie der IR proportional zur Intensität)
    L_ideal = (r_ref ./ distance_ideal).^2;

    % Aktualisiere all_levels für Plot-Grenzen
    all_levels = [];
    for v_idx = 1:numel(results)
        all_levels = [all_levels, results(v_idx).levels];
    end

    fprintf('✓ Alle Werte als relative Energie E/E_ref dargestellt\n');
    fprintf('  (Ideale Kurve: E(r) = E_ref * (r_ref/r)²)\n');
end

%% ---------------- Plot erstellen ----------------
fprintf('\n=== Erstelle Plots ===\n');

% Farben für verschiedene Varianten
colors = lines(numel(varianten));

%% --- PLOT 1: Standard Scatter Plot (ohne Verbindungslinien) ---
fig1 = figure('Position', [100, 100, 1400, 800], 'Visible', 'on');

% Plot für jede Variante
hold on;
legend_entries = {};
plot_handles = [];

for v_idx = 1:numel(results)
    if isempty(results(v_idx).distances)
        continue;
    end

    % Sortiere nach Entfernung für saubere Anzeige
    [sorted_dist, sort_idx] = sort(results(v_idx).distances);
    sorted_levels = results(v_idx).levels(sort_idx);
    sorted_positions = results(v_idx).positions(sort_idx);

    % Plot nur Messpunkte ohne Verbindungslinien
    h = plot(sorted_dist, sorted_levels, 'o', ...
        'Color', colors(v_idx, :), ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', colors(v_idx, :), ...
        'LineWidth', 1.5, ...
        'DisplayName', results(v_idx).variante);
    plot_handles(end+1) = h;
    legend_entries{end+1} = strrep(results(v_idx).variante, '_', ' ');

    % Positionsnummern als Text anzeigen
    for i = 1:length(sorted_dist)
        text(sorted_dist(i), sorted_levels(i), sprintf('  %d', sorted_positions(i)), ...
            'FontSize', 8, 'Color', colors(v_idx, :), ...
            'VerticalAlignment', 'bottom');
    end
end

% Ideale 1/r-Kurve (gestrichelt, schwarz)
h_ideal = plot(distance_ideal, L_ideal, 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal 1/r');
plot_handles(end+1) = h_ideal;
legend_entries{end+1} = sprintf('Ideal 1/r (Halbraum, -6 dB/Abstandsverdopplung)');

% 1/3-Ebene als transparente graue Fläche
y_third_level = max(all_levels) / 3;
x_limits = [0, max_dist * 1.1];
h_third = fill([x_limits(1), x_limits(2), x_limits(2), x_limits(1)], ...
               [y_third_level, y_third_level, y_third_level, y_third_level], ...
               [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
               'DisplayName', '1/3 Pegel-Ebene');
plot_handles(end+1) = h_third;
legend_entries{end+1} = '1/3 Pegel-Ebene';

hold off;
grid on;

% Achsenbeschriftung
xlabel('Entfernung von der Quelle s [m]', 'FontSize', 14, 'FontWeight', 'bold');
if strcmp(darstellung_modus, 'energie')
    % Energie-Darstellung
    if strcmp(analyse_modus, 'summenpegel')
        ylabel('Relative Intensität I/I_{ref}', 'FontSize', 14, 'FontWeight', 'bold');
        title_str = 'Summen-Intensität über Entfernung';
    else
        ylabel(sprintf('Relative Intensität I/I_{ref} bei %.0f Hz', target_frequency), 'FontSize', 14, 'FontWeight', 'bold');
        title_str = sprintf('Intensität (%.0f Hz) über Entfernung', target_frequency);
    end
else
    % Pegel-Darstellung
    if strcmp(analyse_modus, 'summenpegel')
        ylabel('Summenpegel L [dBFS]', 'FontSize', 14, 'FontWeight', 'bold');
        title_str = 'Summenpegel über Entfernung';
    else
        ylabel(sprintf('Terzpegel L [dBFS] bei %.0f Hz', target_frequency), 'FontSize', 14, 'FontWeight', 'bold');
        title_str = sprintf('Terzpegel (%.0f Hz) über Entfernung', target_frequency);
    end
end

title(title_str, 'FontSize', 16, 'FontWeight', 'bold');

% Legende
legend(plot_handles, legend_entries, 'Location', 'northeast', 'FontSize', 11);

% Layout
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'w');

% Achsengrenzen sinnvoll setzen
xlim([0, max_dist * 1.1]);
y_range = max(all_levels) - min(all_levels);
ylim([min(all_levels) - 0.1*y_range, max(all_levels) + 0.1*y_range]);

fprintf('✓ Standard Plot wurde erstellt.\n');

% Plot speichern
if strcmp(darstellung_modus, 'energie')
    if strcmp(analyse_modus, 'summenpegel')
        save_name = 'Summen_Intensitaet_ueber_Entfernung.png';
    else
        save_name = sprintf('Intensitaet_%dHz_ueber_Entfernung.png', target_frequency);
    end
else
    if strcmp(analyse_modus, 'summenpegel')
        save_name = 'Summenpegel_ueber_Entfernung.png';
    else
        save_name = sprintf('Terzpegel_%dHz_ueber_Entfernung.png', target_frequency);
    end
end
if ishandle(fig1) && isvalid(fig1)
    saveas(fig1, save_name);
    fprintf('✓ Plot gespeichert als: %s\n', save_name);
else
    warning('Figure 1 Handle ist ungültig, Plot konnte nicht gespeichert werden.');
end

%% --- PLOT 2: Scatter Plot (2D) für alle Varianten zusammen ---
fig2 = figure('Position', [150, 150, 1400, 800], 'Visible', 'on');
hold on;

% Scatter für jede Variante
for v_idx = 1:numel(results)
    if isempty(results(v_idx).distances)
        continue;
    end

    scatter(results(v_idx).distances, results(v_idx).levels, 150, colors(v_idx, :), ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
        'DisplayName', strrep(results(v_idx).variante, '_', ' '));

    % Positionsnummern anzeigen
    for i = 1:length(results(v_idx).distances)
        text(results(v_idx).distances(i), results(v_idx).levels(i), ...
            sprintf('  %d', results(v_idx).positions(i)), ...
            'FontSize', 9, 'Color', colors(v_idx, :), ...
            'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
    end
end

% Ideale 1/r-Kurve
plot(distance_ideal, L_ideal, 'k--', 'LineWidth', 3, ...
    'DisplayName', 'Ideal 1/r (Halbraum)');

hold off;
grid on;
xlabel('Entfernung von der Quelle s [m]', 'FontSize', 14, 'FontWeight', 'bold');
if strcmp(darstellung_modus, 'energie')
    if strcmp(analyse_modus, 'summenpegel')
        ylabel('Relative Intensität I/I_{ref}', 'FontSize', 14, 'FontWeight', 'bold');
        title_str_scatter = 'Summen-Intensität über Entfernung (Scatter Plot)';
    else
        ylabel(sprintf('Relative Intensität I/I_{ref} bei %.0f Hz', target_frequency), 'FontSize', 14, 'FontWeight', 'bold');
        title_str_scatter = sprintf('Intensität (%.0f Hz) über Entfernung (Scatter Plot)', target_frequency);
    end
else
    if strcmp(analyse_modus, 'summenpegel')
        ylabel('Summenpegel L [dBFS]', 'FontSize', 14, 'FontWeight', 'bold');
        title_str_scatter = 'Summenpegel über Entfernung (Scatter Plot)';
    else
        ylabel(sprintf('Terzpegel L [dBFS] bei %.0f Hz', target_frequency), 'FontSize', 14, 'FontWeight', 'bold');
        title_str_scatter = sprintf('Terzpegel (%.0f Hz) über Entfernung (Scatter Plot)', target_frequency);
    end
end
title(title_str_scatter, 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 11);
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'w');
xlim([0, max_dist * 1.1]);
ylim([min(all_levels) - 0.1*y_range, max(all_levels) + 0.1*y_range]);

% Speichern
if strcmp(darstellung_modus, 'energie')
    if strcmp(analyse_modus, 'summenpegel')
        save_name_scatter = 'Summen_Intensitaet_ueber_Entfernung_Scatter.png';
    else
        save_name_scatter = sprintf('Intensitaet_%dHz_ueber_Entfernung_Scatter.png', target_frequency);
    end
else
    if strcmp(analyse_modus, 'summenpegel')
        save_name_scatter = 'Summenpegel_ueber_Entfernung_Scatter.png';
    else
        save_name_scatter = sprintf('Terzpegel_%dHz_ueber_Entfernung_Scatter.png', target_frequency);
    end
end
if ishandle(fig2) && isvalid(fig2)
    saveas(fig2, save_name_scatter);
    fprintf('✓ Scatter Plot gespeichert als: %s\n', save_name_scatter);
else
    warning('Figure 2 Handle ist ungültig, Plot konnte nicht gespeichert werden.');
end

%% --- PLOT 2b: 2D Plots mit Pfaden für jede Variante separat ---

% Verbindungspfade definieren (gleich wie im 3D-Plot)
line_connections_2d = {
    [13, 14, 15]; ...
    [11, 8]; ...
    [10, 7, 4]; ...
    [6, 3]; ...
    [9, 5, 1]
};

% Erstelle separate 2D-Plots für jede Variante mit Pfaden
for v_idx = 1:numel(results)
    if isempty(results(v_idx).distances)
        continue;
    end

    fig2b = figure('Position', [150 + v_idx*50, 150 + v_idx*50, 1400, 800], 'Visible', 'on');
    hold on;

    % Scatter Plot für diese Variante
    scatter(results(v_idx).distances, results(v_idx).levels, 150, colors(v_idx, :), ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
        'DisplayName', strrep(results(v_idx).variante, '_', ' '));

    % Positionsnummern anzeigen
    for i = 1:length(results(v_idx).distances)
        text(results(v_idx).distances(i), results(v_idx).levels(i), ...
            sprintf('  %d', results(v_idx).positions(i)), ...
            'FontSize', 9, 'Color', colors(v_idx, :), ...
            'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
    end

    % Zeichne Pfade als Regressionskurven
    for line_idx = 1:length(line_connections_2d)
        path_positions = line_connections_2d{line_idx};
        dist_line = [];
        level_line = [];

        % Sammle Punkte entlang des Pfades
        for pos_in_path = path_positions
            pos_result_idx = find(results(v_idx).positions == pos_in_path, 1);
            if ~isempty(pos_result_idx)
                dist_line(end+1) = results(v_idx).distances(pos_result_idx);
                level_line(end+1) = results(v_idx).levels(pos_result_idx);
            end
        end

        % Zeichne Regressionskurve (nur wenn genug Punkte gefunden wurden)
        if length(dist_line) >= 3
            % Sortiere nach Entfernung
            [dist_sorted, sort_idx] = sort(dist_line);
            level_sorted = level_line(sort_idx);

            % Quadratische Regression 2. Ordnung
            p = polyfit(dist_sorted, level_sorted, min(2, length(dist_sorted)-1));

            % Erstelle glatte Kurve
            dist_smooth = linspace(min(dist_sorted), max(dist_sorted), 100);
            level_smooth = polyval(p, dist_smooth);

            % Zeichne Kurve in Variantenfarbe mit Transparenz und gestrichelt
            h_line = plot(dist_smooth, level_smooth, '--', 'Color', colors(v_idx, :), ...
                'LineWidth', 2, 'HandleVisibility', 'off');
            h_line.Color(4) = 0.5; % 50% Deckkraft
        elseif length(dist_line) == 2
            % Bei nur 2 Punkten: einfache gestrichelte Linie
            [dist_sorted, sort_idx] = sort(dist_line);
            level_sorted = level_line(sort_idx);
            h_line = plot(dist_sorted, level_sorted, '--', 'Color', colors(v_idx, :), ...
                'LineWidth', 2, 'HandleVisibility', 'off');
            h_line.Color(4) = 0.5; % 50% Deckkraft
        end
    end

    % Ideale 1/r-Kurve für Halbraum
    plot(distance_ideal, L_ideal, 'k--', 'LineWidth', 3, ...
        'DisplayName', 'Ideal 1/r (Halbraum)');

    hold off;
    grid on;
    xlabel('Entfernung von der Quelle s [m]', 'FontSize', 14, 'FontWeight', 'bold');
    if strcmp(darstellung_modus, 'energie')
        if strcmp(analyse_modus, 'summenpegel')
            ylabel('Relative Intensität I/I_{ref}', 'FontSize', 14, 'FontWeight', 'bold');
            title_str_2d = sprintf('Summen-Intensität über Entfernung - %s', strrep(results(v_idx).variante, '_', ' '));
        else
            ylabel(sprintf('Relative Intensität I/I_{ref} bei %.0f Hz', target_frequency), 'FontSize', 14, 'FontWeight', 'bold');
            title_str_2d = sprintf('Intensität (%.0f Hz) über Entfernung - %s', target_frequency, strrep(results(v_idx).variante, '_', ' '));
        end
    else
        if strcmp(analyse_modus, 'summenpegel')
            ylabel('Summenpegel L [dBFS]', 'FontSize', 14, 'FontWeight', 'bold');
            title_str_2d = sprintf('Summenpegel über Entfernung - %s', strrep(results(v_idx).variante, '_', ' '));
        else
            ylabel(sprintf('Terzpegel L [dBFS] bei %.0f Hz', target_frequency), 'FontSize', 14, 'FontWeight', 'bold');
            title_str_2d = sprintf('Terzpegel (%.0f Hz) über Entfernung - %s', target_frequency, strrep(results(v_idx).variante, '_', ' '));
        end
    end
    title(title_str_2d, 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 11);
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'w');
    xlim([0, max_dist * 1.1]);
    ylim([min(all_levels) - 0.1*y_range, max(all_levels) + 0.1*y_range]);

    % Speichern
    if strcmp(darstellung_modus, 'energie')
        if strcmp(analyse_modus, 'summenpegel')
            save_name_2d_path = sprintf('Summen_Intensitaet_2D_Pfade_%s.png', results(v_idx).variante);
        else
            save_name_2d_path = sprintf('Intensitaet_%dHz_2D_Pfade_%s.png', target_frequency, results(v_idx).variante);
        end
    else
        if strcmp(analyse_modus, 'summenpegel')
            save_name_2d_path = sprintf('Summenpegel_2D_Pfade_%s.png', results(v_idx).variante);
        else
            save_name_2d_path = sprintf('Terzpegel_%dHz_2D_Pfade_%s.png', target_frequency, results(v_idx).variante);
        end
    end
    if ishandle(fig2b) && isvalid(fig2b)
        saveas(fig2b, save_name_2d_path);
        fprintf('✓ 2D Plot mit Pfaden für %s gespeichert als: %s\n', results(v_idx).variante, save_name_2d_path);
    else
        warning('Figure 2b Handle ist ungültig für %s.', results(v_idx).variante);
    end
end

%% --- PLOT 2c: Vergleich beider Varianten mit Differenz ---
if numel(results) == 2
    fig2c = figure('Position', [250, 250, 1400, 900], 'Visible', 'on');

    % Subplot 1: Beide Varianten übereinander mit Differenz-Linien
    subplot(2, 1, 1);
    hold on;

    % Finde gemeinsame Positionen für Differenz-Linien
    common_positions = intersect(results(1).positions, results(2).positions);

    % Zeichne transparente Verbindungslinien zwischen den Varianten für jede Position
    for pos = common_positions
        idx1 = find(results(1).positions == pos, 1);
        idx2 = find(results(2).positions == pos, 1);

        if ~isempty(idx1) && ~isempty(idx2)
            % Vertikale Linie zwischen den beiden Pegeln
            x_pos = results(1).distances(idx1);
            y1 = results(1).levels(idx1);
            y2 = results(2).levels(idx2);

            % Farbe basierend auf Differenz
            if y1 > y2
                line_color = [1 0.7 0.7]; % Rötlich
            else
                line_color = [0.7 0.7 1]; % Bläulich
            end

            plot([x_pos, x_pos], [y1, y2], '-', 'Color', line_color, ...
                'LineWidth', 3, 'HandleVisibility', 'off');
            h = plot([x_pos, x_pos], [y1, y2], '-', 'Color', line_color, ...
                'LineWidth', 3, 'HandleVisibility', 'off');
            h.Color(4) = 0.6; % 60% Transparenz
        end
    end

    for v_idx = 1:2
        if isempty(results(v_idx).distances)
            continue;
        end

        % Scatter Plot
        scatter(results(v_idx).distances, results(v_idx).levels, 150, colors(v_idx, :), ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
            'DisplayName', strrep(results(v_idx).variante, '_', ' '));

        % Positionsnummern
        for i = 1:length(results(v_idx).distances)
            text(results(v_idx).distances(i), results(v_idx).levels(i), ...
                sprintf('  %d', results(v_idx).positions(i)), ...
                'FontSize', 8, 'Color', colors(v_idx, :), ...
                'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
        end
    end

    % Ideale 1/r-Kurve
    plot(distance_ideal, L_ideal, 'k--', 'LineWidth', 2, ...
        'DisplayName', 'Ideal 1/r (Halbraum)');

    hold off;
    grid on;
    xlabel('Entfernung von der Quelle s [m]', 'FontSize', 12, 'FontWeight', 'bold');
    if strcmp(darstellung_modus, 'energie')
        ylabel('Relative Intensität I/I_{ref}', 'FontSize', 12, 'FontWeight', 'bold');
    else
        ylabel('Summenpegel L [dBFS]', 'FontSize', 12, 'FontWeight', 'bold');
    end
    title('Vergleich beider Varianten', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 10);
    set(gca, 'FontSize', 11);
    xlim([0, max_dist * 1.1]);
    ylim([min(all_levels) - 0.1*y_range, max(all_levels) + 0.1*y_range]);

    % Subplot 2: Differenz zwischen den Varianten
    subplot(2, 1, 2);
    hold on;

    % Finde gemeinsame Positionen
    common_positions = intersect(results(1).positions, results(2).positions);

    if ~isempty(common_positions)
        diff_distances = [];
        diff_levels = [];
        diff_positions = [];

        for pos = common_positions
            idx1 = find(results(1).positions == pos, 1);
            idx2 = find(results(2).positions == pos, 1);

            if ~isempty(idx1) && ~isempty(idx2)
                diff_distances(end+1) = results(1).distances(idx1);
                % Differenz: Variante 1 - Variante 2
                diff_levels(end+1) = results(1).levels(idx1) - results(2).levels(idx2);
                diff_positions(end+1) = pos;
            end
        end

        % Plot Differenz
        scatter(diff_distances, diff_levels, 150, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
            'MarkerFaceColor', [0.8 0.3 0.3], ...
            'DisplayName', sprintf('%s - %s', ...
                strrep(results(1).variante, '_', ' '), ...
                strrep(results(2).variante, '_', ' ')));

        % Positionsnummern
        for i = 1:length(diff_distances)
            text(diff_distances(i), diff_levels(i), ...
                sprintf('  %d', diff_positions(i)), ...
                'FontSize', 8, 'Color', [0.8 0.3 0.3], ...
                'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
        end

        % Nulllinie
        plot([0, max_dist * 1.1], [0, 0], 'k-', 'LineWidth', 1.5, ...
            'DisplayName', 'Keine Differenz');

        % Farbige Bereiche für positive/negative Differenzen
        y_max_diff = max(abs(diff_levels)) * 1.2;
        fill([0, max_dist * 1.1, max_dist * 1.1, 0], ...
             [0, 0, y_max_diff, y_max_diff], ...
             [1 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
             'DisplayName', sprintf('%s höher', strrep(results(1).variante, '_', ' ')));
        fill([0, max_dist * 1.1, max_dist * 1.1, 0], ...
             [0, 0, -y_max_diff, -y_max_diff], ...
             [0.9 0.9 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
             'DisplayName', sprintf('%s höher', strrep(results(2).variante, '_', ' ')));

        % Statistik zur Differenz
        mean_diff = mean(diff_levels);
        std_diff = std(diff_levels);
        max_diff = max(abs(diff_levels));

        % Text mit Statistik
        if strcmp(darstellung_modus, 'energie')
            text_str = sprintf('Mittlere Differenz: %.2e\nStandardabweichung: %.2e\nMax. Differenz: %.2e', ...
                mean_diff, std_diff, max_diff);
        else
            text_str = sprintf('Mittlere Differenz: %.2f dB\nStandardabweichung: %.2f dB\nMax. Differenz: %.2f dB', ...
                mean_diff, std_diff, max_diff);
        end
        text(0.02, 0.98, text_str, 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'FontSize', 10, ...
            'BackgroundColor', 'w', 'EdgeColor', 'k');
    end

    hold off;
    grid on;
    xlabel('Entfernung von der Quelle s [m]', 'FontSize', 12, 'FontWeight', 'bold');
    if strcmp(darstellung_modus, 'energie')
        ylabel('Energie-Differenz ΔE/E_{ref}', 'FontSize', 12, 'FontWeight', 'bold');
    else
        ylabel('Pegel-Differenz [dB]', 'FontSize', 12, 'FontWeight', 'bold');
    end
    title(sprintf('Differenz: %s minus %s', ...
        strrep(results(1).variante, '_', ' '), ...
        strrep(results(2).variante, '_', ' ')), ...
        'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 10);
    set(gca, 'FontSize', 11);
    xlim([0, max_dist * 1.1]);

    set(gcf, 'Color', 'w');

    % Speichern
    if strcmp(darstellung_modus, 'energie')
        if strcmp(analyse_modus, 'summenpegel')
            save_name_compare = 'Summen_Intensitaet_Vergleich_mit_Differenz.png';
        else
            save_name_compare = sprintf('Intensitaet_%dHz_Vergleich_mit_Differenz.png', target_frequency);
        end
    else
        if strcmp(analyse_modus, 'summenpegel')
            save_name_compare = 'Summenpegel_Vergleich_mit_Differenz.png';
        else
            save_name_compare = sprintf('Terzpegel_%dHz_Vergleich_mit_Differenz.png', target_frequency);
        end
    end

    if ishandle(fig2c) && isvalid(fig2c)
        saveas(fig2c, save_name_compare);
        fprintf('✓ Vergleichsplot mit Differenz gespeichert als: %s\n', save_name_compare);
    else
        warning('Figure 2c Handle ist ungültig.');
    end
else
    fprintf('ℹ Vergleichsplot übersprungen (benötigt genau 2 Varianten, %d vorhanden)\n', numel(results));
end

%% --- Erstelle 1/r-Ebene Mesh-Daten für 3D Plots ---
% Diese werden für PLOT 2d und PLOT 3 benötigt

x_mesh = linspace(-0.1, max([positions_info.x])*1.1, 30);
y_mesh = linspace(-0.1, max([positions_info.y])*1.1, 30);
[X_mesh, Y_mesh] = meshgrid(x_mesh, y_mesh);
R_mesh = sqrt((X_mesh - source_x).^2 + (Y_mesh - source_y).^2);

% Berechne Werte für die 1/r-Ebene (Halbraum-Ausbreitung)
% An der Quelle (r=0) würde der Wert theoretisch unendlich sein
% Verwende einen physikalisch sinnvollen Mindestabstand basierend auf der kleinsten Messung
min_physical_dist = min(0.05, r_ref * 0.1); % 5cm oder 10% der Referenzentfernung
R_mesh_safe = R_mesh;
R_mesh_safe(R_mesh_safe < min_physical_dist) = min_physical_dist;

if strcmp(darstellung_modus, 'energie')
    % Energie-Modus: E(r) = E_ref * (r_ref/r)^2
    L_mesh = (r_ref ./ R_mesh_safe).^2;
    L_source = (r_ref / min_physical_dist)^2;
    fprintf('1/r²-Ebene (Halbraum, Energie): Mindestabstand = %.3f m, E_source/E_ref = %.2f\n', min_physical_dist, L_source);
else
    % Pegel-Modus: L(r) = L_ref - 20*log10(r/r_ref)
    L_mesh = L_ref - 20*log10(R_mesh_safe / r_ref);
    L_source = L_ref - 20*log10(min_physical_dist / r_ref);
    fprintf('1/r-Ebene (Halbraum, Pegel): Mindestabstand = %.3f m, L_source = %.1f dBFS\n', min_physical_dist, L_source);
end

%% --- PLOT 2d: 3D Vergleichsplot beider Varianten ---
if numel(results) == 2
    fig2d = figure('Position', [300, 300, 1400, 900], 'Visible', 'on');
    hold on;

    % Finde gemeinsame Positionen
    common_positions = intersect(results(1).positions, results(2).positions);

    % Sammle x,y-Koordinaten für beide Varianten
    for v_idx = 1:2
        if isempty(results(v_idx).distances)
            continue;
        end

        x_coords = [];
        y_coords = [];
        for p_idx = 1:length(results(v_idx).positions)
            pos = results(v_idx).positions(p_idx);
            pos_info_idx = find([positions_info.pos] == pos);
            if ~isempty(pos_info_idx)
                x_coords(end+1) = positions_info(pos_info_idx).x;
                y_coords(end+1) = positions_info(pos_info_idx).y;
            end
        end

        % 3D Scatter Plot für diese Variante
        scatter3(x_coords, y_coords, results(v_idx).levels, 150, colors(v_idx, :), ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
            'DisplayName', strrep(results(v_idx).variante, '_', ' '));

        % Positionsnummern in 3D anzeigen
        for i = 1:length(x_coords)
            text(x_coords(i), y_coords(i), results(v_idx).levels(i), ...
                sprintf('  P%d', results(v_idx).positions(i)), ...
                'FontSize', 9, 'Color', colors(v_idx, :), ...
                'FontWeight', 'bold');
        end
    end

    % Zeichne transparente Verbindungslinien zwischen den Varianten
    for pos = common_positions
        idx1 = find(results(1).positions == pos, 1);
        idx2 = find(results(2).positions == pos, 1);

        if ~isempty(idx1) && ~isempty(idx2)
            % Finde x,y Koordinaten
            pos_info_idx = find([positions_info.pos] == pos);
            if ~isempty(pos_info_idx)
                x_pos = positions_info(pos_info_idx).x;
                y_pos = positions_info(pos_info_idx).y;
                z1 = results(1).levels(idx1);
                z2 = results(2).levels(idx2);

                % Vertikale Linie zwischen den beiden Pegeln
                if z1 > z2
                    line_color = [1 0.5 0.5]; % Rot: Variante 1 höher
                else
                    line_color = [0.5 0.5 1]; % Blau: Variante 2 höher
                end

                h = plot3([x_pos, x_pos], [y_pos, y_pos], [z1, z2], '-', ...
                    'Color', line_color, 'LineWidth', 4, 'HandleVisibility', 'off');
                h.Color(4) = 0.5; % 50% Transparenz
            end
        end
    end

    % Quellenposition als Kugel/Kreis markieren (auf der 1/r-Ebene)
    scatter3(source_x, source_y, L_source, 300, 'r', 'o', ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2.5, ...
        'DisplayName', 'Quelle Q1');
    text(source_x, source_y, L_source, '  Q1', ...
        'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');

    % Zeichne die 1/r-Ebene als transparentes Mesh (Halbraum)
    surf(X_mesh, Y_mesh, L_mesh, 'FaceAlpha', 0.2, 'EdgeColor', 'k', ...
        'EdgeAlpha', 0.05, 'FaceColor', [0.7 0.7 0.7], ...
        'DisplayName', 'Ideal 1/r-Ebene (Halbraum)');

    hold off;
    grid on;
    xlabel('x-Position [m]', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('y-Position [m]', 'FontSize', 13, 'FontWeight', 'bold');
    if strcmp(darstellung_modus, 'energie')
        if strcmp(analyse_modus, 'summenpegel')
            zlabel('Relative Intensität E/E_{ref}', 'FontSize', 13, 'FontWeight', 'bold');
            title_str_3d_comp = 'Vergleich beider Varianten im Raum - Energie (3D)';
        else
            zlabel(sprintf('Relative Intensität E/E_{ref} bei %.0f Hz', target_frequency), 'FontSize', 13, 'FontWeight', 'bold');
            title_str_3d_comp = sprintf('Vergleich beider Varianten (%.0f Hz) im Raum - Energie (3D)', target_frequency);
        end
    else
        if strcmp(analyse_modus, 'summenpegel')
            zlabel('Summenpegel L [dBFS]', 'FontSize', 13, 'FontWeight', 'bold');
            title_str_3d_comp = 'Vergleich beider Varianten im Raum (3D)';
        else
            zlabel(sprintf('Terzpegel L [dBFS] bei %.0f Hz', target_frequency), 'FontSize', 13, 'FontWeight', 'bold');
            title_str_3d_comp = sprintf('Vergleich beider Varianten (%.0f Hz) im Raum (3D)', target_frequency);
        end
    end
    title(title_str_3d_comp, 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 11);
    view(45, 30);
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'w');

    % Speichern
    if strcmp(darstellung_modus, 'energie')
        if strcmp(analyse_modus, 'summenpegel')
            save_name_3d_comp = 'Summen_Intensitaet_3D_Vergleich_beider_Varianten.png';
        else
            save_name_3d_comp = sprintf('Intensitaet_%dHz_3D_Vergleich_beider_Varianten.png', target_frequency);
        end
    else
        if strcmp(analyse_modus, 'summenpegel')
            save_name_3d_comp = 'Summenpegel_3D_Vergleich_beider_Varianten.png';
        else
            save_name_3d_comp = sprintf('Terzpegel_%dHz_3D_Vergleich_beider_Varianten.png', target_frequency);
        end
    end

    if ishandle(fig2d) && isvalid(fig2d)
        saveas(fig2d, save_name_3d_comp);
        fprintf('✓ 3D Vergleichsplot gespeichert als: %s\n', save_name_3d_comp);
    else
        warning('Figure 2d Handle ist ungültig.');
    end
else
    fprintf('ℹ 3D Vergleichsplot übersprungen (benötigt genau 2 Varianten, %d vorhanden)\n', numel(results));
end

%% --- PLOT 3: Scatter3 Plot für jede Variante separat (3D mit x, y, Pegel) ---

% Verbindungslinien von Quelle zu Messpunkten als Regressionskurven
% Q-P13-14-15, Q-11-8, Q-P10-P7-P4, Q-P6-P3, Q-P9-P5-P1
line_connections = {
    [13, 14, 15]; ...
    [11, 8]; ...
    [10, 7, 4]; ...
    [6, 3]; ...
    [9, 5, 1]
};

% Erstelle separate 3D-Plots für jede Variante
for v_idx = 1:numel(results)
    if isempty(results(v_idx).distances)
        continue;
    end

    fig3 = figure('Position', [200 + v_idx*50, 200 + v_idx*50, 1400, 900], 'Visible', 'on');
    hold on;

    % Hole x,y-Koordinaten für die Messpositionen dieser Variante
    x_coords = [];
    y_coords = [];
    for p_idx = 1:length(results(v_idx).positions)
        pos = results(v_idx).positions(p_idx);
        pos_info_idx = find([positions_info.pos] == pos);
        if ~isempty(pos_info_idx)
            x_coords(end+1) = positions_info(pos_info_idx).x;
            y_coords(end+1) = positions_info(pos_info_idx).y;
        end
    end

    % 3D Scatter Plot für diese Variante
    scatter3(x_coords, y_coords, results(v_idx).levels, 150, colors(v_idx, :), ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
        'DisplayName', strrep(results(v_idx).variante, '_', ' '));

    % Positionsnummern in 3D anzeigen
    for i = 1:length(x_coords)
        text(x_coords(i), y_coords(i), results(v_idx).levels(i), ...
            sprintf('  P%d', results(v_idx).positions(i)), ...
            'FontSize', 9, 'Color', colors(v_idx, :), ...
            'FontWeight', 'bold');
    end

    % Quellenposition als Kugel/Kreis markieren (auf der 1/r-Ebene)
    scatter3(source_x, source_y, L_source, 300, 'r', 'o', ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2.5, ...
        'DisplayName', 'Quelle Q1');
    text(source_x, source_y, L_source, '  Q1', ...
        'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');

    % Zeichne Regressionskurven für diese Variante
    for line_idx = 1:length(line_connections)
        path_positions = line_connections{line_idx};
        x_line = source_x;
        y_line = source_y;
        z_line = L_source;

        for pos_in_path = path_positions
            pos_result_idx = find(results(v_idx).positions == pos_in_path, 1);
            if ~isempty(pos_result_idx)
                pos_info_idx = find([positions_info.pos] == pos_in_path);
                if ~isempty(pos_info_idx)
                    x_line(end+1) = positions_info(pos_info_idx).x;
                    y_line(end+1) = positions_info(pos_info_idx).y;
                    z_line(end+1) = results(v_idx).levels(pos_result_idx);
                end
            end
        end

        % Zeichne Regressionskurve (nur wenn genug Punkte gefunden wurden)
        if length(x_line) >= 3
            dist_along_path = zeros(size(x_line));
            for i = 2:length(x_line)
                dist_along_path(i) = dist_along_path(i-1) + ...
                    sqrt((x_line(i)-x_line(i-1))^2 + (y_line(i)-y_line(i-1))^2);
            end

            p_x = polyfit(dist_along_path, x_line, min(2, length(x_line)-1));
            p_y = polyfit(dist_along_path, y_line, min(2, length(y_line)-1));
            p_z = polyfit(dist_along_path, z_line, min(2, length(z_line)-1));

            dist_smooth = linspace(0, dist_along_path(end), 100);
            x_smooth = polyval(p_x, dist_smooth);
            y_smooth = polyval(p_y, dist_smooth);
            z_smooth = polyval(p_z, dist_smooth);

            h_line = plot3(x_smooth, y_smooth, z_smooth, '--', 'Color', colors(v_idx, :), ...
                'LineWidth', 2, 'HandleVisibility', 'off');
            h_line.Color(4) = 0.5; % 50% Deckkraft
        elseif length(x_line) == 2
            h_line = plot3(x_line, y_line, z_line, '--', 'Color', colors(v_idx, :), ...
                'LineWidth', 2, 'HandleVisibility', 'off');
            h_line.Color(4) = 0.5; % 50% Deckkraft
        end
    end

    % Zeichne die 1/r-Ebene als transparentes Mesh (Halbraum)
    surf(X_mesh, Y_mesh, L_mesh, 'FaceAlpha', 0.3, 'EdgeColor', 'k', ...
        'EdgeAlpha', 0.1, 'FaceColor', [0.7 0.7 0.7], ...
        'DisplayName', 'Ideal 1/r-Ebene (Halbraum)');

    hold off;
    grid on;
    xlabel('x-Position [m]', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('y-Position [m]', 'FontSize', 13, 'FontWeight', 'bold');
    if strcmp(darstellung_modus, 'energie')
        if strcmp(analyse_modus, 'summenpegel')
            zlabel('Relative Intensität E/E_{ref}', 'FontSize', 13, 'FontWeight', 'bold');
            title_str_3d = sprintf('Summen-Intensität im Raum - %s', strrep(results(v_idx).variante, '_', ' '));
        else
            zlabel(sprintf('Relative Intensität E/E_{ref} bei %.0f Hz', target_frequency), 'FontSize', 13, 'FontWeight', 'bold');
            title_str_3d = sprintf('Intensität (%.0f Hz) im Raum - %s', target_frequency, strrep(results(v_idx).variante, '_', ' '));
        end
    else
        if strcmp(analyse_modus, 'summenpegel')
            zlabel('Summenpegel L [dBFS]', 'FontSize', 13, 'FontWeight', 'bold');
            title_str_3d = sprintf('Summenpegel im Raum - %s', strrep(results(v_idx).variante, '_', ' '));
        else
            zlabel(sprintf('Terzpegel L [dBFS] bei %.0f Hz', target_frequency), 'FontSize', 13, 'FontWeight', 'bold');
            title_str_3d = sprintf('Terzpegel (%.0f Hz) im Raum - %s', target_frequency, strrep(results(v_idx).variante, '_', ' '));
        end
    end
    title(title_str_3d, 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 11);
    view(45, 30);
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'w');

    % Speichern
    if strcmp(darstellung_modus, 'energie')
        if strcmp(analyse_modus, 'summenpegel')
            save_name_3d = sprintf('Summen_Intensitaet_3D_%s.png', results(v_idx).variante);
        else
            save_name_3d = sprintf('Intensitaet_%dHz_3D_%s.png', target_frequency, results(v_idx).variante);
        end
    else
        if strcmp(analyse_modus, 'summenpegel')
            save_name_3d = sprintf('Summenpegel_3D_%s.png', results(v_idx).variante);
        else
            save_name_3d = sprintf('Terzpegel_%dHz_3D_%s.png', target_frequency, results(v_idx).variante);
        end
    end
    if ishandle(fig3) && isvalid(fig3)
        saveas(fig3, save_name_3d);
        fprintf('✓ 3D Plot für %s gespeichert als: %s\n', results(v_idx).variante, save_name_3d);
    else
        warning('Figure 3 Handle ist ungültig für %s.', results(v_idx).variante);
    end
end

%% ---------------- Statistische Auswertung ----------------
fprintf('\n=== Statistische Auswertung ===\n');
fprintf('=================================================================\n');

% Sammle alle Abweichungen für Gesamtstatistik
all_deviations = [];
all_std_values = [];

for v_idx = 1:numel(results)
    if isempty(results(v_idx).distances)
        continue;
    end

    fprintf('\n--- %s ---\n', strrep(results(v_idx).variante, '_', ' '));

    % Berechne Abweichung von idealer 1/r-Kurve
    L_ideal_at_meas = L_ref - 20*log10(results(v_idx).distances / r_ref);
    deviation = results(v_idx).levels - L_ideal_at_meas;

    % Speichere für Gesamtstatistik
    all_deviations = [all_deviations, deviation];
    std_dev = std(deviation);
    all_std_values(end+1) = std_dev;

    % Statistische Kennwerte
    fprintf('  Anzahl Messpunkte: %d\n', length(deviation));
    fprintf('  \n');
    fprintf('  ABWEICHUNG VON IDEALER 1/r-DÄMPFUNG:\n');
    fprintf('  • Mittlere Abweichung:        %+.2f dB\n', mean(deviation));
    fprintf('  • Standardabweichung (STD):    %.2f dB  ← Hauptkennwert!\n', std_dev);
    fprintf('  • Varianz:                     %.2f dB²\n', var(deviation));
    fprintf('  • RMS-Abweichung:              %.2f dB\n', sqrt(mean(deviation.^2)));
    fprintf('  • Max. positive Abweichung:   %+.2f dB (weniger Dämpfung als Freifeld)\n', max(deviation));
    fprintf('  • Max. negative Abweichung:   %+.2f dB (mehr Dämpfung als Freifeld)\n', min(deviation));
    fprintf('  • Spannweite (Range):          %.2f dB\n', max(deviation) - min(deviation));
    fprintf('  \n');

    % Berechne effektiven Dämpfungsexponenten durch Regression
    % L = a - b*log10(r)  => b sollte 20 für ideale 1/r sein
    if length(results(v_idx).distances) > 2
        log_dist = log10(results(v_idx).distances);
        p = polyfit(log_dist, results(v_idx).levels, 1);
        effective_exponent = -p(1);
        deviation_from_ideal_exp = effective_exponent - 20;

        fprintf('  DÄMPFUNGSEXPONENT:\n');
        fprintf('  • Gemessener Exponent:         %.2f dB/Dekade\n', effective_exponent);
        fprintf('  • Idealer Exponent (1/r):      20.00 dB/Dekade\n');
        fprintf('  • Abweichung:                 %+.2f dB/Dekade (%.1f%%)\n', ...
            deviation_from_ideal_exp, (deviation_from_ideal_exp/20)*100);
        fprintf('  • Gemessen (pro Verdopplung):  %.2f dB\n', effective_exponent * log10(2));
        fprintf('  • Ideal (pro Verdopplung):     6.02 dB\n');

        % Bestimmtheitsmaß R² der Regression
        y_fit = polyval(p, log_dist);
        SS_res = sum((results(v_idx).levels - y_fit).^2);
        SS_tot = sum((results(v_idx).levels - mean(results(v_idx).levels)).^2);
        R_squared = 1 - SS_res/SS_tot;
        fprintf('  • R² (Güte der Anpassung):     %.4f\n', R_squared);
    end

    fprintf('  \n');

    % Speichere Abweichungen in results für spätere Verwendung
    results(v_idx).deviations = deviation;
    results(v_idx).std_deviation = std_dev;
end

% Gesamtstatistik über alle Varianten
if numel(results) > 1 && ~isempty(all_deviations)
    fprintf('\n=================================================================\n');
    fprintf('--- GESAMTSTATISTIK ÜBER ALLE VARIANTEN ---\n');
    fprintf('  Gesamtanzahl Messpunkte: %d\n', length(all_deviations));
    fprintf('  \n');
    fprintf('  • Mittlere STD über Varianten: %.2f dB\n', mean(all_std_values));
    fprintf('  • Globale STD (alle Punkte):   %.2f dB\n', std(all_deviations));
    fprintf('  • Globaler Mittelwert:        %+.2f dB\n', mean(all_deviations));
    fprintf('  • Globale RMS-Abweichung:      %.2f dB\n', sqrt(mean(all_deviations.^2)));
    fprintf('=================================================================\n');
end

fprintf('\n✓ Statistische Analyse abgeschlossen.\n');
fprintf('\nINTERPRETATION:\n');
fprintf('• STD nahe 0 dB  → Sehr gute Übereinstimmung mit Halbraum-Ausbreitung\n');
fprintf('• STD < 2 dB     → Gute Übereinstimmung, geringe Reflexionen\n');
fprintf('• STD > 3 dB     → Deutliche Abweichungen, Raumeinflüsse sichtbar\n');
fprintf('• Positive Mittelwerte → Weniger Dämpfung als erwartet (Reflexionen)\n');
fprintf('• Negative Mittelwerte → Mehr Dämpfung als erwartet (Absorption)\n');


%% ========================================================================
%%                        HILFSFUNKTIONEN
%% ========================================================================

%% Berechnet Summenpegel für eine Position
function L_summe = getSummenpegel(variante, position, dataDir, allFiles, fs, f_terz, FS_global)
    % Holt Terzpegel für alle Bänder
    terzpegel = getTerzpegel(variante, position, dataDir, allFiles, fs, f_terz, FS_global);

    if all(isnan(terzpegel))
        L_summe = NaN;
        return;
    end

    % Energetische Summation über alle Frequenzbänder
    % L_gesamt = 10*log10(sum(10^(L_i/10)))
    linear_values = 10.^(terzpegel / 10);
    L_summe = 10 * log10(sum(linear_values, 'omitnan'));
end

%% Berechnet Terzpegel für alle Bänder einer Position
function L_dBFS = getTerzpegel(variante, position, dataDir, allFiles, fs, f_terz, FS_global)

    L_dBFS = NaN(1, length(f_terz)); % Initialisieren mit NaN

    % Passende Datei finden
    posStr = num2str(position);
    fname = '';

    % Suche nach verschiedenen Namensmustern
    % Standard: Variante_X,Pos_Y.mat oder Variante_X_neu,Pos_Y.mat
    patterns = {
        sprintf('%s,Pos_%s.mat', variante, posStr),
        sprintf('%s_neu,Pos_%s.mat', variante, posStr),
        sprintf('%s,Pos%s.mat', variante, posStr)
    };

    for p = 1:length(patterns)
        for i = 1:numel(allFiles)
            if strcmp(allFiles{i}, patterns{p})
                fname = allFiles{i};
                break;
            end
        end
        if ~isempty(fname), break; end
    end

    % Flexiblere Suche mit RegEx falls nicht gefunden
    if isempty(fname)
        pattern = ['^' variante '.*Pos[_,]?' posStr '([^\d]|\.mat$)'];
        for i = 1:numel(allFiles)
             if ~isempty(regexp(allFiles{i}, pattern, 'once'))
                fname = allFiles{i};
                break;
            end
        end
    end

    if isempty(fname)
        return;
    end

    % Lade und verarbeite Datei
    try
        S = load(fullfile(dataDir, fname));
        ir = extractIR(S);

        if isempty(ir)
            warning('Datei %s enthält keinen erkennbaren IR-Vektor.', fname);
            return;
        end

        % Lundeby-Truncation
        [ir_trunc, ~, ~, ~, ~, ~] = truncateIR(ir);
        N = length(ir_trunc);

        % FFT berechnen
        IR_fft = fft(ir_trunc, N);
        freq = (0:N-1) * (fs / N);

        % Terzband-Energie berechnen
        f_lower = f_terz / 2^(1/6);
        f_upper = f_terz * 2^(1/6);
        IR_fft_abs2 = abs(IR_fft).^2;

        for k = 1:length(f_terz)
            if f_upper(k) >= fs/2, continue; end
            idx = (freq >= f_lower(k)) & (freq <= f_upper(k));
            if ~any(idx), continue; end

            total_energy = sum(IR_fft_abs2(idx)) / N;
            rms_band = sqrt(total_energy);
            L_dBFS(k) = 20*log10((rms_band + eps) / FS_global);
        end

    catch ME
        warning('Fehler bei der Verarbeitung von %s: %s', fname, ME.message);
        return;
    end
end

%% Extrahiert Impulsantwort aus .mat Struktur
function ir = extractIR(S)
    ir = [];
    if isfield(S,'RiR') && ~isempty(S.RiR)
        ir = double(S.RiR(:));
    elseif isfield(S,'RIR') && ~isempty(S.RIR)
        ir = double(S.RIR(:));
    else
        fns = fieldnames(S);
        for f = 1:numel(fns)
            fname = fns{f};
            if startsWith(fname, '__'), continue; end
            v = S.(fname);
            if isnumeric(v) && numel(v) > 1000
                ir = double(v(:));
                return;
            end
        end
    end
    if ~isempty(ir) && numel(ir) < 2
        ir = [];
    end
end

%% Lundeby-Truncation zur Rauschunterdrückung
function [ir_trunc, start_idx, end_idx, E_ratio, SNR_dB, dynamic_range_dB] = truncateIR(ir)
    N_original = length(ir);
    ir_abs = abs(ir);
    max_amp = max(ir_abs);

    % Rauschpegel aus letzten 10% schätzen
    noise_samples = ir_abs(end-round(N_original*0.1):end);
    noise_level = mean(noise_samples);
    noise_rms = std(noise_samples);

    % Schwellwert für Signalende
    threshold = max((noise_level + 3*noise_rms) * 10, max_amp * 0.001);
    sig_idx = find(ir_abs > threshold, 1, 'last');

    if isempty(sig_idx) || sig_idx < 100
        end_idx = N_original;
    else
        end_idx = min(N_original, round(sig_idx * 1.2));
    end

    % Signalbeginn finden (5% des Maximums)
    start_threshold = max_amp * 0.05;
    start_idx = find(ir_abs > start_threshold, 1, 'first');
    if isempty(start_idx) || start_idx < 1
        start_idx = 1;
    end

    ir_trunc = ir(start_idx:end_idx);

    % Qualitätsmetriken
    E_original = sum(ir.^2);
    E_truncated = sum(ir_trunc.^2);
    E_ratio = E_truncated / E_original * 100;

    signal_rms = sqrt(mean(ir_trunc.^2));
    SNR_dB = 20*log10(signal_rms / (noise_rms + eps));

    dynamic_range_dB = 20*log10(max_amp / (noise_level + eps));
end

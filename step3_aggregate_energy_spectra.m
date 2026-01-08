%% step3_aggregate_energy_spectra.m
%  Aggregiert Summenenergie und Terzpegel aller Varianten
%  Erstellt eine Excel-Datei mit Übersicht und Spektren.

clear; clc;

% Konfiguration
procDir = 'processed';
outFile = 'Zusammenfassung_Energie_Terz.xlsx';

if ~exist(procDir, 'dir')
    error('Ordner "processed" nicht gefunden. Bitte erst step1 ausführen.');
end

% Dateien finden
files = dir(fullfile(procDir, 'Proc_*.mat'));
if isempty(files)
    error('Keine verarbeiteten Dateien gefunden.');
end

fprintf('Verarbeite %d Dateien...\n', length(files));

% Container
data_rows = {};
spectra_matrix = [];
freq_vec = [];

for i = 1:length(files)
    % Laden
    tmp = load(fullfile(files(i).folder, files(i).name), 'Result');
    R = tmp.Result;
    
    % Metadaten
    varName = R.meta.variante;
    posName = R.meta.position;
    if isfield(R.meta, 'type') && strcmp(R.meta.type, 'Source')
        posName = 'Quelle';
    end
    
    % Werte
    E_lin = R.time.metrics.energy;      % Lineare Energie (sum x^2)
    L_sum = R.freq.sum_level;           % Summenpegel dB
    L_terz = R.freq.terz_dbfs;          % Terzspektrum
    
    if isempty(freq_vec)
        freq_vec = R.freq.f_center;
    end
    
    % Zeile hinzufügen
    data_rows(end+1, :) = {files(i).name, varName, posName, E_lin, L_sum};
    spectra_matrix(end+1, :) = L_terz;
end

% Tabelle 1: Übersicht
T_main = cell2table(data_rows, 'VariableNames', {'Datei', 'Variante', 'Position', 'Energie_Linear', 'Summenpegel_dB'});

% Sortieren und im Command Window ausgeben
T_main = sortrows(T_main, {'Variante', 'Position'});
disp('--- Übersicht: Summenpegel der einzelnen Positionen ---');
disp(T_main(:, {'Variante', 'Position', 'Summenpegel_dB'}));

% Tabelle 2: Spektren
% Spaltennamen generieren (z.B. Hz_50, Hz_63...)
colNames = arrayfun(@(f) sprintf('Hz_%d', round(f)), freq_vec, 'UniformOutput', false);
colNames = strrep(colNames, '.', '_'); % Sicherheitshalber
T_spec = array2table(spectra_matrix, 'VariableNames', colNames);

% Alles in eine Tabelle
T_complete = [T_main, T_spec];

% Speichern
writetable(T_complete, fullfile(procDir, outFile), 'Sheet', 'Alle_Daten');
fprintf('Daten gespeichert in: %s (Sheet: Alle_Daten)\n', fullfile(procDir, outFile));

% --- Aggregation pro Variante (Summe über alle Positionen) ---
variants = unique(T_main.Variante);
agg_rows = {};

for k = 1:length(variants)
    v = variants{k};
    % Filter: Nur diese Variante und KEINE Quelle (nur Empfängerpositionen)
    mask = strcmp(T_main.Variante, v) & ~strcmp(T_main.Position, 'Quelle');
    
    if any(mask)
        % Summe der linearen Energien
        total_energy = sum(T_main.Energie_Linear(mask));
        
        % Durchschnittlicher Summenpegel (arithmetisch in dB)
        avg_level = mean(T_main.Summenpegel_dB(mask));
        
        % Durchschnittliches Spektrum
        avg_spec = mean(spectra_matrix(mask, :), 1);
        
        agg_rows(end+1, :) = [{v, total_energy, avg_level}, num2cell(avg_spec)];
    end
end

if ~isempty(agg_rows)
    T_agg = cell2table(agg_rows, 'VariableNames', [{'Variante', 'Gesamt_Energie_Linear', 'Durchschnitt_Pegel_dB'}, colNames]);
    writetable(T_agg, fullfile(procDir, outFile), 'Sheet', 'Varianten_Vergleich');
    fprintf('Aggregierte Daten gespeichert in Sheet: Varianten_Vergleich\n');
end

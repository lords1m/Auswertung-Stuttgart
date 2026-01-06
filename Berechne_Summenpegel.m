%% ============================================================ 
%  Berechnung und Speicherung der Summenpegel
% 
%  Dieses Skript analysiert alle .mat-Dateien im 'data'-Ordner,
%  berechnet für jede Messung (definiert durch Variante und Position)
%  den Summenpegel aus den Terzspektren (dBFS) und speichert
%  die Ergebnisse in einer Excel-Tabelle.
% 
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

%% ---------------- Einstellungen ---------------- 

dataDir = 'data'; % Datenordner mit allen .mat-Dateien
outputDir = 'Excel'; % Ausgabeordner für die Excel-Datei
outputFileName = 'Summenpegel_Aller_Varianten.xlsx';

% --- Physikalische Parameter (aus Terzpegel_DBFs_einzeln.m übernommen) --- 
fs = 500e3; % 500 kHz Abtastrate

% Normgerechte 1/3-Oktav-Mittenfrequenzen (IEC 61260)
f_terz = double([ ...
    560 625 710 800 900 1000 1120 1250 1400 1600 1800 2000 ...
    2240 2500 2800 3150 3550 4000 4500 5000 5600 6300 7100 8000 ...
    9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 ...
    28000 31500 35500 40000 45000 50000 56000 63000 71000 80000 ...
    90000 100000 112000 126000 ]);

%% ---------------- Setup & Dateiprüfung ---------------- 
if ~exist(dataDir, 'dir'), error('Datenordner "%s" nicht gefunden!', dataDir); end
if ~exist(outputDir,'dir'), mkdir(outputDir); end

dirInfo = dir(fullfile(dataDir, '*.mat'));
matFiles = {dirInfo.name};
if isempty(matFiles), error('Keine .mat-Dateien im Ordner "%s" gefunden!', dataDir); end

%% ---------------- Globale Referenz (FS) für dBFS finden ---------------- 
% Notwendig, damit alle Spektren auf derselben Skala vergleichbar sind
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
if FS_global == 0, error('Keine gültigen Impulsantworten gefunden. FS_Global konnte nicht bestimmt werden.'); end
fprintf('✓ Globale Referenz (FS) ermittelt: %g\n', FS_global);


%% ---------------- Summenpegel für jede Datei berechnen ---------------- 
fprintf('\nBerechne Summenpegel für alle Messungen...\n');

% Initialisiere eine Tabelle für die Ergebnisse
resultsTable = table('Size', [numel(matFiles), 4], ...
                     'VariableTypes', {'string', 'string', 'string', 'double'}, ...
                     'VariableNames', {'Datei', 'Variante', 'Position', 'Summenpegel_dBFS'});

for i = 1:numel(matFiles)
    fname = matFiles{i};
    fprintf('  Verarbeite: %s\n', fname);
    
    % Extrahiere Variante und Position aus dem Dateinamen
    [variante, posStr] = parseFileName(fname);
    
    % Terzpegel berechnen (modifizierte Logik aus getTerzpegel)
    L_dBFS = calculateTerzpegelForFile(fullfile(dataDir, fname), fs, f_terz, FS_global);
    
    if all(isnan(L_dBFS))
        warning('Keine gültigen Terzpegel für %s berechnet.', fname);
        summenpegel = NaN;
    else
        % Summenpegel aus Terzpegeln berechnen
        summenpegel = 10 * log10(sum(10.^(L_dBFS(~isnan(L_dBFS)) / 10)));
        fprintf('    -> Summenpegel: %.2f dBFS\n', summenpegel);
    end
    
    % Ergebnisse in Tabelle speichern
    resultsTable(i, :) = {fname, variante, posStr, summenpegel};
end

%% ---------------- Ergebnisse in Excel speichern ---------------- 
% Sortiere die Tabelle für bessere Lesbarkeit
resultsTable = sortrows(resultsTable, {'Variante', 'Position'});

outputFilePath = fullfile(outputDir, outputFileName);
try
    writetable(resultsTable, outputFilePath, 'Sheet', 'Summenpegel');
    fprintf('\nErgebnisse erfolgreich gespeichert in:\n  %s\n', outputFilePath);
catch ME
    error('Fehler beim Speichern der Excel-Datei: %s', ME.message);
end


%% ======================================================================== 
%  HILFSFUNKTIONEN
%  ======================================================================== 

function L_dBFS = calculateTerzpegelForFile(filePath, fs, f_terz, FS_global)
    L_dBFS = NaN(1, length(f_terz)); % Initialisieren mit NaN
    
    try
        S = load(filePath);
        ir = extractIR(S);
        
        if isempty(ir)
            warning('Datei %s enthält keinen erkennbaren IR-Vektor.', filePath);
            return;
        end
        
        [ir_trunc, ~, ~, ~, ~, ~] = truncateIR(ir);
        N = length(ir_trunc);

        IR_fft = fft(ir_trunc, N);
        freq = (0:N-1) * (fs / N);
        
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
        warning('Fehler bei der Verarbeitung von %s: %s', filePath, ME.message);
        return;
    end
end

function [variante, posStr] = parseFileName(fname)
    % Extrahiert Variante und Position aus dem Dateinamen
    % Beispiel: 'Variante_1_neu,Pos_10.mat' -> 'Variante_1_neu', '10'
    
    variante = "N/A";
    posStr = "N/A";
    
    % RegEx, um Variante und Position zu erfassen
    % Updated pattern to be more robust
    pattern = '^(.*?),Pos_([A-Za-z0-9_]+)\.mat$';
    tokens = regexp(fname, pattern, 'tokens');
    
    if ~isempty(tokens)
        variante = string(tokens{1}{1});
        posStr = string(tokens{1}{2});
    else
        warning('Konnte Variante/Position aus Dateinamen nicht extrahieren: %s', fname);
    end
end

% Die folgenden Hilfsfunktionen sind direkt aus Terzpegel_DBFs_einzeln.m kopiert
function ir = extractIR(S)
    ir = [];
    if isfield(S,'RiR') && ~isempty(S.RiR), ir = double(S.RiR(:));
    elseif isfield(S,'RIR') && ~isempty(S.RIR), ir = double(S.RIR(:));
    else
        fns = fieldnames(S);
        for f = 1:numel(fns)
            fname = fns{f};
            if startsWith(fname, '__'), continue; end
            v = S.(fname);
            if isnumeric(v) && numel(v) > 1000, ir = double(v(:)); return; end
        end
    end
    if ~isempty(ir) && numel(ir) < 2, ir = []; end
end

function [ir_trunc, start_idx, end_idx, E_ratio, SNR_dB, dynamic_range_dB] = truncateIR(ir)
    N_original = length(ir);
    ir_abs = abs(ir);
    max_amp = max(ir_abs);
    if max_amp == 0, ir_trunc = ir; start_idx=1; end_idx=N_original; E_ratio=100; SNR_dB=0; dynamic_range_dB=0; return; end

    noise_samples = ir_abs(end-round(N_original*0.1):end);
    noise_level = mean(noise_samples);
    noise_rms = std(noise_samples);

    threshold = max((noise_level + 3*noise_rms) * 10, max_amp * 0.001);
    sig_idx = find(ir_abs > threshold, 1, 'last');

    if isempty(sig_idx) || sig_idx < 100, end_idx = N_original;
    else, end_idx = min(N_original, round(sig_idx * 1.2)); end

    start_threshold = max_amp * 0.05;
    start_idx = find(ir_abs > start_threshold, 1, 'first');
    if isempty(start_idx) || start_idx < 1, start_idx = 1; end

    ir_trunc = ir(start_idx:end_idx);
    
    E_original = sum(ir.^2);
    E_truncated = sum(ir_trunc.^2);
    if E_original == 0, E_ratio = 100; else, E_ratio = E_truncated / E_original * 100; end

    signal_rms = sqrt(mean(ir_trunc.^2));
    SNR_dB = 20*log10(signal_rms / (noise_rms + eps));
    dynamic_range_dB = 20*log10(max_amp / (noise_level + eps));
end

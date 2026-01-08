%% step1_process_data.m
%  Zentrales Skript: Liest Raw-Data -> Berechnet -> Speichert Processed-Data
%  UPDATE: Optimiert für Formate wie 'Variante_2_Pos_15.mat'
clear; clc; close all;

% Pfade setzen und Helper laden
scriptDir = fileparts(mfilename('fullpath'));
if ~isempty(scriptDir), cd(scriptDir); end
addpath('functions'); 

% --- Config ---
dataDir = 'data';
procDir = 'processed';
fs = 500e3; % Abtastrate

if ~exist(dataDir, 'dir'), error('Ordner "data" fehlt!'); end
if ~exist(procDir, 'dir'), mkdir(procDir); end

% Ordner für separate Funktions-Ausgaben erstellen
dirTime = fullfile(procDir, 'Time_Domain');
dirFreq = fullfile(procDir, 'Frequency_Domain');
if ~exist(dirTime, 'dir'), mkdir(dirTime); end
if ~exist(dirFreq, 'dir'), mkdir(dirFreq); end

% 1. Globale Referenz finden
files = dir(fullfile(dataDir, '*.mat'));
fprintf('Schritt 1: Scanne %d Dateien für globale Referenz...\n', length(files));

FS_global = 0;
for i = 1:length(files)
    try
        filepath = fullfile(files(i).folder, files(i).name);
        % Kurzer Load für Max-Check
        S = load(filepath);
        ir = extract_ir(S);
        if ~isempty(ir)
             FS_global = max(FS_global, max(abs(ir)));
        end
    catch
    end
end
if FS_global == 0, FS_global = 1; end
fprintf('Globale Referenz (Max Amp): %.5f\n', FS_global);

% 2. Verarbeitungsschleife
fprintf('Schritt 2: Verarbeitung und Speicherung...\n');
summary_data = {}; 

for i = 1:length(files)
    filepath = fullfile(files(i).folder, files(i).name);
    
    % A) Parsen (liefert z.B. Variante='Variante_2', Position='15')
    [S, meta] = load_and_parse_file(filepath);
    
    % Überspringen bei Fehlern
    if isempty(S) || isempty(meta.variante)
        fprintf('  [SKIP] %s\n', files(i).name);
        continue;
    end
    
    % B) IR extrahieren
    ir_raw = extract_ir(S);
    if isempty(ir_raw), continue; end
    
    % C) Truncation
    [ir_trunc, metrics] = truncate_ir(ir_raw);
    
    % Namen generieren für separate Files
    if strcmp(meta.type, 'Source')
        nameTag = sprintf('%s_Quelle', meta.variante);
    else
        nameTag = sprintf('%s_Pos%s', meta.variante, meta.position);
    end
    
    % Speichern: Truncated IR & Metrics (Time Domain)
    save(fullfile(dirTime, ['Time_' nameTag '.mat']), 'ir_trunc', 'metrics', 'meta');
    
    % D) Spektrum & Summenpegel
    [L_terz, L_sum, f_center] = calc_terz_spectrum(ir_trunc, fs, FS_global);
    
    % Speichern: Spektrum (Frequency Domain)
    save(fullfile(dirFreq, ['Spec_' nameTag '.mat']), 'L_terz', 'L_sum', 'f_center', 'meta');
    
    % E) Ergebnis-Struktur
    Result = struct();
    Result.meta = meta;
    Result.meta.fs = fs;
    Result.meta.FS_global_used = FS_global;
    
    Result.time.ir = ir_trunc;
    Result.time.metrics = metrics;
    
    Result.freq.f_center = f_center;
    Result.freq.terz_dbfs = L_terz;
    Result.freq.sum_level = L_sum;
    
    % F) Speichern
    % Erzeugt standardisierte Namen: 'Proc_Variante_2_Pos15.mat'
    if strcmp(meta.type, 'Source')
        saveName = sprintf('Proc_%s_Quelle.mat', meta.variante);
    else
        % Hier wird Pos%s genutzt, um auch 'Z1' zu unterstützen
        saveName = sprintf('Proc_%s_Pos%s.mat', meta.variante, meta.position);
    end
    
    save(fullfile(procDir, saveName), 'Result');
    
    % Eintrag für Tabelle
    summary_data(end+1,:) = {meta.variante, meta.position, L_sum, metrics.snr_db, saveName};
    
    if mod(i, 5) == 0, fprintf('.'); end
end

% 3. Zusammenfassung speichern
if ~isempty(summary_data)
    summary_table = cell2table(summary_data, ...
        'VariableNames', {'Variante', 'Position', 'SumLevel', 'SNR', 'File'});
    
    summary_table = sortrows(summary_table, {'Variante', 'Position'});
    
    save(fullfile(procDir, 'Summary_Database.mat'), 'summary_table');
    writetable(summary_table, fullfile(procDir, 'Summary.xlsx'));
    
    % --- NEU: Durchschnittsspektren (Average) ---
    fprintf('\nSchritt 3: Erstelle energetisch gemittelte Spektren (Average)...\n');
    unique_vars = unique(summary_table.Variante);
    
    for i = 1:length(unique_vars)
        v = unique_vars{i};
        % Nur Positionen (keine Quelle)
        mask = strcmp(summary_table.Variante, v) & ~strcmp(summary_table.Position, 'Quelle');
        subset = summary_table(mask, :);
        
        if isempty(subset), continue; end
        
        % Template laden (für Metadaten)
        first_file = fullfile(procDir, subset.File{1});
        tmp = load(first_file);
        R_tmpl = tmp.Result;
        
        % Akkumulatoren
        sum_E_terz = 0;
        sum_E_total = 0;
        n = height(subset);
        
        for k = 1:n
            D = load(fullfile(procDir, subset.File{k}));
            % Terzpegel -> Lineare Energie (relativ zu FS)
            E_terz = 10.^(D.Result.freq.terz_dbfs / 10);
            E_terz(~isfinite(E_terz)) = 0;
            sum_E_terz = sum_E_terz + E_terz;
            
            % Summenpegel -> Lineare Energie
            E_tot = 10^(D.Result.freq.sum_level / 10);
            if ~isfinite(E_tot), E_tot = 0; end
            sum_E_total = sum_E_total + E_tot;
        end
        
        % Mittelung
        avg_E_terz = sum_E_terz / n;
        avg_E_total = sum_E_total / n;
        
        % Result bauen
        Result = R_tmpl;
        Result.meta.position = 'Average';
        Result.meta.type = 'Average';
        Result.time.ir = zeros(100,1); % Dummy IR, da spektral gemittelt
        Result.time.metrics.energy = avg_E_total * (Result.meta.FS_global_used^2); % Abs. Energie
        Result.time.metrics.snr_db = NaN;
        Result.time.metrics.idx_start = 1;
        Result.time.metrics.idx_end = 1;
        Result.time.metrics.energy_total = NaN;
        Result.time.metrics.energy_share = NaN;
        
        Result.freq.terz_dbfs = 10 * log10(avg_E_terz + eps);
        Result.freq.sum_level = 10 * log10(avg_E_total + eps);
        
        saveName = sprintf('Proc_%s_Average.mat', v);
        save(fullfile(procDir, saveName), 'Result');
        fprintf('  -> %s erstellt.\n', saveName);
    end
    
    fprintf('\nFertig! %d Dateien verarbeitet.\nDaten in "%s".\n', height(summary_table), procDir);
else
    warning('Keine Daten verarbeitet.');
end

% --- Lokale Helper Funktionen ---
function [S, meta] = load_and_parse_file(filepath)
    [~, fname, ~] = fileparts(filepath);
    S = load(filepath);
    meta = struct();
    meta.filename = fname;
    
    % Regex für Variante und Position
    tokens = regexp(fname, '^(.*?)[_,]Pos[_,]?(\w+)', 'tokens', 'once', 'ignorecase');
    if ~isempty(tokens)
        meta.variante = tokens{1};
        meta.position = tokens{2};
        meta.type = 'Receiver';
    else
        tokens = regexp(fname, '^(.*?)[_,]Quelle', 'tokens', 'once', 'ignorecase');
        if ~isempty(tokens)
            meta.variante = tokens{1};
            meta.position = 'Q1';
            meta.type = 'Source';
        else
            meta.variante = 'Unknown';
            meta.position = '0';
            meta.type = 'Unknown';
        end
    end
end

function ir = extract_ir(S)
    ir = [];
    if isfield(S,'RiR') && ~isempty(S.RiR), ir = double(S.RiR(:));
    elseif isfield(S,'RIR') && ~isempty(S.RIR), ir = double(S.RIR(:));
    elseif isfield(S,'aufn') && ~isempty(S.aufn), ir = double(S.aufn(:));
    else
        fns = fieldnames(S);
        for f = 1:numel(fns)
            fname = fns{f};
            if startsWith(fname, '__'), continue; end
            v = S.(fname);
            if isnumeric(v) && numel(v) > 1000, ir = double(v(:)); return; end
        end
    end
end
%% ============================================================
%  Terzpegel-Auswertung für ausgewählte Einzelmessungen
%  
%  Dieses Skript ermöglicht es, gezielt einzelne oder mehrere
%  Messpunkte (Position und Variante) auszuwählen und deren
%  Terzspektren gemeinsam in einem einzigen Plot darzustellen.
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

%% ---------------- Einstellungen (anpassen) ----------------

% --- Messungen zum Plotten auswählen ---
% Fügen Sie hier jede gewünschte Kombination hinzu.
% - Für eine Einzelmessung: { 'Varianten-Name', Positions-Nummer }
% - Für eine ganze Variante: { 'Varianten-Name', [] } (leere Klammern)
messungen = { ...
    {'Variante_1', []}, ...       % Plot-Beispiel: Alle Positionen von Variante 1
    {'Variante_2', 1}, ...      % Plot-Beispiel: Nur Position 1 von Variante 2
};

% --- Allgemeine Einstellungen ---
dataDir = 'data'; % Datenordner mit allen .mat-Dateien
outputPlotDir = 'Plots'; % Ausgabeordner für den Plot
selectedPositions = 1:14; % Wenn eine ganze Variante geplottet wird, werden diese Positionen verwendet

% --- Physikalische & Plot-Parameter ---
fs = 500e3; % 500 kHz Abtastrate
xtick_vals = [4000 5000 10000 20000 50000 60000];
xtick_labels = {'4k','5k','10k','20k','50k','60k'};
f_lim = [4000 60000];

% Normgerechte 1/3-Oktav-Mittenfrequenzen (IEC 61260)
f_terz = double([ ...
    560 625 710 800 900 1000 1120 1250 1400 1600 1800 2000 ...
    2240 2500 2800 3150 3550 4000 4500 5000 5600 6300 7100 8000 ...
    9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 ...
    28000 31500 35500 40000 45000 50000 56000 63000 71000 80000 ...
    90000 100000 112000 126000 ]);

%% ---------------- Setup & Dateiprüfung ----------------
if ~exist(dataDir, 'dir'), error('Datenordner "%s" nicht gefunden!', dataDir); end
if ~exist(outputPlotDir,'dir'), mkdir(outputPlotDir); end

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

%% ---------------- Terzpegel berechnen & plotten ----------------
all_L_dBFS = []; % Sammel-Array für die Y-Achsen-Limits
legend_entries = {};

fig = figure('Visible','on','Position',[100,100,1200,600]);
hold on;

fprintf('\nVerarbeite ausgewählte Messungen...\n');
for i = 1:numel(messungen)
    variante = messungen{i}{1};
    position_spec = messungen{i}{2};
    
    positions_to_plot = position_spec;
    if isempty(position_spec)
        % Fall: Ganze Variante plotten
        positions_to_plot = selectedPositions;
        fprintf('  Verarbeite ganze Variante: %s (Positionen %d-%d)\n', variante, min(positions_to_plot), max(positions_to_plot));
    end
    
    for j = 1:numel(positions_to_plot)
        pos = positions_to_plot(j);
        
        if isempty(position_spec) % Nur loggen, wenn es eine Einzelmessung ist
             % Kein Log-Output für jede einzelne Position bei "ganzer Variante"
        else
             fprintf('  Lese ein: %s, Pos %d\n', variante, pos);
        end

        % Terzpegel für diese eine Messung berechnen
        L_dBFS = getTerzpegel(variante, pos, dataDir, matFiles, fs, f_terz, FS_global);

        if ~all(isnan(L_dBFS))
            stairs(f_terz, L_dBFS, 'LineWidth', 1.5);
            all_L_dBFS = [all_L_dBFS; L_dBFS]; % Für Skalierung der Y-Achse
            legend_entries{end+1} = sprintf('%s, Pos %d', variante, pos);
        else
            warning('Keine gültigen Daten für %s, Pos %d erhalten. Wird im Plot ausgelassen.', variante, pos);
        end
    end
end

if isempty(all_L_dBFS)
    close(fig);
    error('Keine einzige der ausgewählten Messungen konnte verarbeitet werden. Es wird kein Plot erstellt.');
end

hold off;
grid on;
set(gca, 'XScale', 'log');

% Achsen und Titel
title('Vergleich ausgewählter Terzpegel', 'FontSize', 14);
xlabel('Frequenz [Hz]', 'FontSize', 12);
ylabel('Pegel [dBFS]', 'FontSize', 12);

% Achsen-Limits
xlim(f_lim);
xticks(xtick_vals);
xticklabels(xtick_labels);
validData = all_L_dBFS(~isnan(all_L_dBFS));
if ~isempty(validData)
    y_range = [floor(min(validData)/10)*10, ceil(max(validData)/10)*10];
    ylim(y_range);
end

% Legende
legend(legend_entries, 'Location', 'southwest', 'Interpreter', 'none');
set(gcf, 'Color', 'w');

% Speichern
outputFileName = fullfile(outputPlotDir, 'Vergleich_Terzpegel_Einzeln.png');
saveas(fig, outputFileName);
saveas(fig, strrep(outputFileName, '.png', '.fig'));

fprintf('\nPlot wurde erfolgreich erstellt und gespeichert:\n  %s\n', outputFileName);


%% ========================================================================
%  HILFSFUNKTIONEN
%  ========================================================================

%% ---------------- Hilfsfunktion zur Terzpegelberechnung (Einzelmessung) ----------------
function L_dBFS = getTerzpegel(variante, position, dataDir, allFiles, fs, f_terz, FS_global)
    
    L_dBFS = NaN(1, length(f_terz)); % Initialisieren mit NaN

    % Passende Datei finden
    posStr = num2str(position);
    fname = '';
    
    % Suche nach exaktem Match zuerst (effizienter)
    searchPattern_exact = sprintf('%s_Pos%s.mat', variante, posStr);
    for i = 1:numel(allFiles)
        if strcmp(allFiles{i}, searchPattern_exact)
            fname = allFiles{i};
            break;
        end
    end
    
    % Wenn nicht gefunden, flexiblere Suche mit RegEx
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
        return; % Keine Warnung hier, wird in der aufrufenden Funktion behandelt
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

%% ---------------- Hilfsfunktionen aus anderen Skripten ----------------
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
%% ============================================================
%  Terzpegel-Auswertung (1/3 Oktav) aus RIR
%  Ergebnis: dBFS-Werte pro Terzband
%  Erweiterungen:
%   - automatische Referenzsuche (global über ausgewählte Varianten, _alt ausgeschlossen)
%   - selektive Auswahl von Varianten und Positionen
%   - Ausgabe: alle Plots in Plots/, alle Excels in Excel/
%   - Plot-Modus: 'absolute' oder 'difference' (Differenz gegen positionsmittleren Pegel)
%   - Alle .mat-Dateien in einem gemeinsamen 'data'-Ordner
% ============================================================

clear;
clc;

% Stelle sicher, dass wir im richtigen Verzeichnis sind
scriptDir = fileparts(mfilename('fullpath'));
if ~isempty(scriptDir)
    cd(scriptDir);
end
fprintf('Arbeitsverzeichnis: %s\n', pwd);

%% ---------------- Einstellungen (anpassen) ----------------
% Datenordner mit allen .mat-Dateien
dataDir = 'data';

% Varianten zum Verarbeiten (leer = automatisch alle außer '_alt' suchen)
selectedVariants = { 'Variante_1'};  % z.B. {'Variante_1', 'Variante_2'} oder {} für Auto
excludePattern = '_alt';  % Ausschlussmuster für Varianten

% Positionen (Messpunkte) auswählen
selectedPositions = 1:14;   % z.B. [1 3 5] oder 1:14

% Ausgabeordner
outputPlotDir = 'Plots';
outputExcelDir = 'Excel';

% Sampling
fs = 500e3;   % 500 kHz Abtastrate

% Normgerechte 1/3-Oktav-Mitttenfrequenzen (IEC 61260)
f_terz = double([ ...
    560 625 710 800 900 1000 1120 1250 1400 1600 1800 2000 ...
    2240 2500 2800 3150 3550 4000 4500 5000 5600 6300 7100 8000 ...
    9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 ...
    28000 31500 35500 40000 45000 50000 56000 63000 71000 80000 ...
    90000 100000 112000 126000 ]);

%% ---------------- Prüfe data-Ordner ----------------
if ~exist(dataDir, 'dir')
    error('Datenordner "%s" nicht gefunden!', dataDir);
end

%% ---------------- Dateien einlesen und Varianten ermitteln ----------------
% Alle .mat-Dateien im data-Ordner finden (mit dir - effizienter als Java)
dirInfo = dir(fullfile(dataDir, '*.mat'));
matFiles = {dirInfo.name};

if isempty(matFiles)
    error('Keine .mat-Dateien im Ordner "%s" gefunden!', dataDir);
end

fprintf('Gefundene .mat-Dateien: %d\n', numel(matFiles));

% Varianten aus Dateinamen extrahieren (mit unique statt Map)
variantNames = cell(numel(matFiles), 1);
validIdx = false(numel(matFiles), 1);

for i = 1:numel(matFiles)
    fname = matFiles{i};

    % Extrahiere Variantenname vor dem ersten Komma oder 'Pos'
    tokens = regexp(fname, '^([^,]+?)(?:_neu)?[,_]Pos', 'tokens', 'once');
    if ~isempty(tokens)
        variantName = tokens{1};

        % Filtern nach excludePattern
        if isempty(excludePattern) || ~contains(variantName, excludePattern)
            variantNames{i} = variantName;
            validIdx(i) = true;
        end
    end
end

variantNames = unique(variantNames(validIdx));

% Wenn user spezifische Varianten gesetzt hat -> filter auf diese Liste
if ~isempty(selectedVariants)
    variantNames = intersect(variantNames, selectedVariants, 'stable');
end

if isempty(variantNames)
    error('Keine Varianten zum Verarbeiten gefunden (nach Filter).');
end

fprintf('\nVerarbeite Varianten:\n');
fprintf('  - %s\n', variantNames{:});

% Vorberechnung: Terzbandgrenzen
nFreq = length(f_terz);
f_lower = f_terz / 2^(1/6);
f_upper = f_terz * 2^(1/6);

% Erstelle File-Lookup-Map für schnelleren Zugriff
fileMap = buildFileMap(matFiles, variantNames, selectedPositions);

%% ---------------- Globale Referenz (FS) über alle ausgewählten Dateien finden ----------------
FS_global = 0;
checkedFiles = 0;

fprintf('\nSuche globale Referenz in %d Varianten...\n', numel(variantNames));

for vi = 1:numel(variantNames)
    variantName = variantNames{vi};
    fprintf('  Prüfe Variante: %s\n', variantName);

    for pos = selectedPositions
        fname = getFileFromMap(fileMap, variantName, pos);
        if isempty(fname)
            continue;
        end

        filename = fullfile(dataDir, fname);

        % Lade und verarbeite Datei
        try
            S = load(filename);
            ir = extractIR(S);

            if ~isempty(ir)
                % Lundeby-Truncation auch für globale Referenz
                [ir, ~, ~, ~, ~, ~] = truncateIR(ir);
                N = numel(ir);
                H_mag = abs(fft(ir, N));
                FS_global = max(FS_global, max(H_mag(1:floor(N/2)+1)));
                checkedFiles = checkedFiles + 1;
            end
        catch
            continue;
        end
    end
end

if FS_global == 0 || checkedFiles == 0
    error('Keine gültigen Impulsantworten gefunden. FS_Global = 0');
end
fprintf('✓ Globale Referenz (FS) ermittelt: %g (Spektral-Maximum) aus %d Dateien\n', FS_global, checkedFiles);

%% ---------------- Ensure output dirs ----------------
if ~exist(outputPlotDir,'dir'), mkdir(outputPlotDir); end
if ~exist(outputExcelDir,'dir'), mkdir(outputExcelDir); end

%% ---------------- Verarbeitung pro Variante ----------------
for vi = 1:numel(variantNames)
    variantName = variantNames{vi};
    fprintf('\n========================================\n');
    fprintf('Verarbeite %s...\n', variantName);
    fprintf('========================================\n');

    % Speichergrößen: Positionen dynamisch
    nPos = numel(selectedPositions);

    L_dBFS = NaN(nPos, nFreq);
    H_all = cell(nPos,1);
    freq_all = cell(nPos,1);

    for pi = 1:nPos
        pos = selectedPositions(pi);

        % Suche passende Datei
        fname = getFileFromMap(fileMap, variantName, pos);
        if isempty(fname)
            warning('Keine Datei für Pos %d in %s gefunden. Position übersprungen.', pos, variantName);
            continue;
        end

        filename = fullfile(dataDir, fname);
        S = load(filename);
        ir = extractIR(S);

        if isempty(ir)
            warning('Datei %s enthält keinen erkennbaren IR-Vektor. Übersprungen.', filename);
            continue;
        end

        % --- Lundeby-Truncation: Finde Rauschgrenze ---
        ir_original = ir;
        N_original = length(ir);
        [ir, start_idx, end_idx, E_ratio, SNR_dB, dynamic_range_dB] = truncateIR(ir);
        N = length(ir);

        % Logging
        fprintf('  Position %02d: Start=%d, Ende=%d (von %d), Länge=%d, Energie=%.2f%%, SNR=%.1f dB, DR=%.1f dB\n', ...
                pos, start_idx, end_idx, N_original, N, E_ratio, SNR_dB, dynamic_range_dB);

        % Warnung und Plot bei geringer Energie
        if E_ratio < 60
            fprintf('  ⚠ WARNUNG: Position %02d hat nur %.2f%% Energie - Plot wird erstellt!\n', pos, E_ratio);
            plotLowEnergyIR(ir_original, ir, start_idx, end_idx, variantName, pos, E_ratio, fs, SNR_dB, dynamic_range_dB);
        end

        IR_fft = fft(ir);
        freq = (0:N-1) * (fs / N);

        nHalf = floor(N/2)+1;
        H_mag = abs(IR_fft(1:nHalf));
        freq_pos = freq(1:nHalf);

        % dBFS relativ zur globalen Referenz
        H_dBFS = 20*log10((H_mag + eps) / FS_global);

        % Terzband-Auswertung (vektorisiert wo möglich)
        IR_fft_abs2 = abs(IR_fft).^2;

        for k = 1:nFreq
            if f_upper(k) >= fs/2
                continue;  % L_dBFS bleibt NaN
            end

            idx = (freq >= f_lower(k)) & (freq <= f_upper(k));
            if ~any(idx)
                continue;  % L_dBFS bleibt NaN
            end

            total_energy = sum(IR_fft_abs2(idx)) / N;
            rms_band = sqrt(total_energy);
            L_dBFS(pi,k) = 20*log10((rms_band + eps) / FS_global);
        end

        H_all{pi} = H_dBFS;
        freq_all{pi} = freq_pos;
    end

    % Mittelwert pro Terzband (über die ausgewählten Positionen)
    valid_mask = ~isnan(L_dBFS);
    if any(valid_mask(:))
        L_mean_dBFS = 10*log10(sum(10.^(L_dBFS/10) .* valid_mask, 1) ./ max(1, sum(valid_mask,1)));
    else
        L_mean_dBFS = NaN(1,nFreq);
    end

    % Summen-Terzpegel (energetische Addition über alle Positionen)
    if any(valid_mask(:))
        L_sum_dBFS = 10*log10(sum(10.^(L_dBFS/10) .* valid_mask, 1));
    else
        L_sum_dBFS = NaN(1,nFreq);
    end

    %% ---------------- Excel-Export ----------------
    rowNames = compose('Pos_%02d', selectedPositions);
    colNames = compose('F%.0f', f_terz);

    T_LdBFS = array2table(L_dBFS, 'RowNames', rowNames, 'VariableNames', colNames);

    excelFile = fullfile(outputExcelDir, sprintf('%s_Terzpegel_dBFS.xlsx', variantName));
    try
        writetable(T_LdBFS, excelFile, 'WriteRowNames', true);
        fprintf('Excel geschrieben: %s\n', excelFile);
    catch ME
        warning(ME.identifier, 'Konnte Excel nicht schreiben: %s', ME.message);
    end

    %% ---------------- Plots ----------------
    % Y-Achsen-Grenzen (global für diese Variante)
    validData = L_dBFS(~isnan(L_dBFS));
    if ~isempty(validData)
        y_range_terz = [floor(min(validData)/10)*10, ceil(max(validData)/10)*10];
    else
        y_range_terz = [-100 0];
    end

    % Plot-Einstellungen (konstant für alle Plots)
    xtick_vals = [4000 5000 10000 20000 50000 60000];
    xtick_labels = {'4k','5k','10k','20k','50k','60k'};
    f_lim = [4000 60000];

    % Terzpegel-Plots
    for pi = 1:nPos
        pos = selectedPositions(pi);
        if all(isnan(L_dBFS(pi,:))), continue; end

        fig = figure('Visible','off','Position',[100,100,1000,500]);
        stairs(f_terz, L_dBFS(pi,:), 'LineWidth',2, 'Color',[0 0.4470 0.7410]);
        grid on; set(gca,'XScale','log');
        xlabel('Frequenz [Hz]'); ylabel('Pegel [dBFS]');
        title(sprintf('%s - Position %02d', variantName, pos));
        xlim(f_lim);
        ylim(y_range_terz);
        xticks(xtick_vals);
        xticklabels(xtick_labels);

        filename = fullfile(outputPlotDir, sprintf('Terzpegel_%s_Pos_%02d.png', variantName, pos));
        saveas(fig, filename);
        saveas(fig, strrep(filename,'.png','.fig'));
        close(fig);
    end

    % Übertragungsfunktionen
    validH = cellfun(@(x) ~isempty(x), H_all);
    if any(validH)
        % Finde Min/Max über alle gültigen Übertragungsfunktionen
        % (ohne cell2mat, da unterschiedliche Längen möglich)
        y_min_H = inf;
        y_max_H = -inf;
        for h_idx = find(validH)'
            y_min_H = min(y_min_H, min(H_all{h_idx}));
            y_max_H = max(y_max_H, max(H_all{h_idx}));
        end
        y_range_H = [floor(y_min_H/10)*10, ceil(y_max_H/10)*10];
    else
        y_range_H = [-120 0];
    end

    for pi = 1:nPos
        pos = selectedPositions(pi);
        if isempty(H_all{pi}), continue; end

        fig = figure('Visible','off','Position',[100,100,800,500]);
        semilogx(freq_all{pi}, H_all{pi}, '-', 'LineWidth',1.5, 'Color',[0.85 0.33 0.10]);
        grid on;
        xlabel('Frequenz [Hz]'); ylabel('Pegel [dBFS]');
        title(sprintf('Uebertragungsfunktion - %s - Pos %02d', variantName, pos));
        xlim(f_lim);
        ylim(y_range_H);
        xticks(xtick_vals);
        xticklabels(xtick_labels);

        filename = fullfile(outputPlotDir, sprintf('Uebertragungsfunktion_%s_Pos_%02d%s.png', variantName, pos, modeSuffix));
        saveas(fig, filename);
        saveas(fig, strrep(filename,'.png','.fig'));
        close(fig);
    end

    % Summen-Terzpegel Plot
    fig = figure('Visible','off','Position',[100,100,1000,500]);
    stairs(f_terz, L_sum_dBFS, 'LineWidth',2, 'Color',[0.85 0.33 0.10]);
    grid on; set(gca,'XScale','log');
    xlabel('Frequenz [Hz]'); ylabel('Summen-Pegel [dB]');
    title(sprintf('%s - Summen-Terzpegel (alle Positionen)', variantName));
    xlim(f_lim);
    xticks(xtick_vals);
    xticklabels(xtick_labels);

    % Y-Achse basierend auf Summen-Daten
    validSum = L_sum_dBFS(~isnan(L_sum_dBFS));
    if ~isempty(validSum)
        ylim([floor(min(validSum)/10)*10, ceil(max(validSum)/10)*10]);
    end

    filename = fullfile(outputPlotDir, sprintf('Terzpegel_%s_Summe%s.png', variantName, modeSuffix));
    saveas(fig, filename);
    saveas(fig, strrep(filename,'.png','.fig'));
    close(fig);

    % Mittelwert-Spektrum (Übertragungsfunktion gemittelt über alle Positionen)
    if any(validH)
        % Finde gemeinsame Frequenz-Auflösung (nehme die mit den meisten Punkten)
        max_len = 0;
        ref_idx = 0;
        for pi = 1:nPos
            if ~isempty(H_all{pi}) && length(H_all{pi}) > max_len
                max_len = length(H_all{pi});
                ref_idx = pi;
            end
        end

        if ref_idx > 0
            % Verwende Frequenz-Vektor der Referenz-Übertragungsfunktion
            freq_mean = freq_all{ref_idx};
            H_mean = zeros(size(freq_mean));
            count = 0;

            % Interpoliere alle Übertragungsfunktionen auf gemeinsame Frequenzen und mittele
            for pi = 1:nPos
                if ~isempty(H_all{pi})
                    % Interpoliere auf gemeinsame Frequenzen
                    H_interp = interp1(freq_all{pi}, H_all{pi}, freq_mean, 'linear', NaN);
                    % Addiere gültige Werte
                    valid_idx = ~isnan(H_interp);
                    H_mean(valid_idx) = H_mean(valid_idx) + H_interp(valid_idx);
                    count = count + 1;
                end
            end

            % Mittelwert berechnen
            if count > 0
                H_mean = H_mean / count;

                % Plot erstellen
                fig = figure('Visible','off','Position',[100,100,1000,500]);
                semilogx(freq_mean, H_mean, '-', 'LineWidth',2, 'Color',[0.4660 0.6740 0.1880]);
                grid on;
                xlabel('Frequenz [Hz]'); ylabel('Pegel [dBFS]');
                title(sprintf('%s - Mittelwert-Spektrum (alle Positionen)', variantName));
                xlim(f_lim);
                ylim(y_range_H);
                xticks(xtick_vals);
                xticklabels(xtick_labels);

                filename = fullfile(outputPlotDir, sprintf('Spektrum_%s_Mittelwert%s.png', variantName, modeSuffix));
                saveas(fig, filename);
                saveas(fig, strrep(filename,'.png','.fig'));
                close(fig);

                fprintf('  → Mittelwert-Spektrum erstellt: %s\n', filename);
            end
        end
    end

    fprintf('%s abgeschlossen: Excel + Plots gespeichert (inkl. Summen-Terzpegel + Mittelwert-Spektrum).\n', variantName);
end

disp('========================================');
disp('Alle ausgewählten Varianten erfolgreich verarbeitet!');
disp('========================================');

%% ---------------- Hilfsfunktionen ----------------
function fileMap = buildFileMap(matFiles, variantNames, positions)
    % Erstellt eine Map für schnellen Dateizugriff
    % fileMap{variant_idx, pos_idx} = filename

    nVariants = numel(variantNames);
    nPos = numel(positions);
    fileMap = cell(nVariants, nPos);

    for i = 1:numel(matFiles)
        fname = matFiles{i};

        % Finde zugehörige Variante
        for vi = 1:nVariants
            if ~startsWith(fname, variantNames{vi})
                continue;
            end

            % Finde Position
            for pi = 1:nPos
                posStr = num2str(positions(pi));
                pattern = ['Pos[_,]' posStr '([^\d]|\.mat$)'];
                if ~isempty(regexp(fname, pattern, 'once'))
                    fileMap{vi, pi} = fname;
                    break;
                end
            end
            break;
        end
    end
end

function fname = getFileFromMap(fileMap, variantName, pos)
    % Holt Dateinamen aus der Map
    fname = '';

    [nVariants, nPos] = size(fileMap);

    % Finde Varianten-Index
    vi = find(strcmp(variantName, fileMap(1:nVariants, 1)), 1);
    if isempty(vi)
        % Lineare Suche als Fallback
        for v = 1:nVariants
            for p = 1:nPos
                if ~isempty(fileMap{v, p}) && startsWith(fileMap{v, p}, variantName)
                    vi = v;
                    break;
                end
            end
            if ~isempty(vi), break; end
        end
    end

    if isempty(vi), return; end

    % Finde Position-Index (lineare Suche durch Positionen)
    for pi = 1:nPos
        if ~isempty(fileMap{vi, pi})
            posStr = num2str(pos);
            pattern = ['Pos[_,]' posStr '([^\d]|\.mat$)'];
            if ~isempty(regexp(fileMap{vi, pi}, pattern, 'once'))
                fname = fileMap{vi, pi};
                return;
            end
        end
    end
end

function ir = extractIR(S)
    % Extrahiert die Impulsantwort aus der geladenen Struktur

    ir = [];

    % Versuche bekannte Feldnamen
    if isfield(S,'RiR') && ~isempty(S.RiR)
        ir = double(S.RiR(:));
    elseif isfield(S,'RIR') && ~isempty(S.RIR)
        ir = double(S.RIR(:));
    else
        % Fallback: erstes numerisches Feld mit vielen Elementen
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

    % Validierung
    if ~isempty(ir) && numel(ir) < 2
        ir = [];
    end
end

function [ir_trunc, start_idx, end_idx, E_ratio, SNR_dB, dynamic_range_dB] = truncateIR(ir)
    % Lundeby-Truncation: Findet Nutz-Signal-Bereich
    % Gibt zurück: truncierte IR, Start-Index, End-Index, Energie-Verhältnis (%), SNR (dB), Dynamikbereich (dB)

    N_original = length(ir);
    ir_abs = abs(ir);
    max_amp = max(ir_abs);

    % Rauschpegel aus letzten 10% schätzen
    noise_samples = ir_abs(end-round(N_original*0.1):end);
    noise_level = mean(noise_samples);
    noise_rms = std(noise_samples);

    % Finde letzten signifikanten Peak (mindestens 10x über Rauschpegel)
    threshold = max((noise_level + 3*noise_rms) * 10, max_amp * 0.001);  % mindestens -60 dB
    sig_idx = find(ir_abs > threshold, 1, 'last');

    if isempty(sig_idx) || sig_idx < 100
        end_idx = N_original;  % Falls keine Truncation möglich
    else
        % Sicherheitsmarge: 20% zusätzlich
        end_idx = min(N_original, round(sig_idx * 1.2));
    end

    % Finde Startzeitpunkt (erster signifikanter Peak)
    start_threshold = max_amp * 0.05;  % 5% des Maximums
    start_idx = find(ir_abs > start_threshold, 1, 'first');
    if isempty(start_idx) || start_idx < 1
        start_idx = 1;
    end

    % Truncierte Impulsantwort
    ir_trunc = ir(start_idx:end_idx);

    % Energie-Verhältnis
    E_original = sum(ir.^2);
    E_truncated = sum(ir_trunc.^2);
    E_ratio = E_truncated / E_original * 100;

    % SNR Berechnung: Signal RMS vs Noise RMS
    signal_rms = sqrt(mean(ir_trunc.^2));
    SNR_dB = 20*log10(signal_rms / (noise_rms + eps));

    % Nutzbare Dynamik: Peak zu Noise Floor
    dynamic_range_dB = 20*log10(max_amp / (noise_level + eps));
end

function plotLowEnergyIR(ir_original, ir_trunc, start_idx, end_idx, variantName, pos, E_ratio, fs, SNR_dB, dynamic_range_dB)
    % Erstellt einen Plot für IRs mit geringer Energie zur Diagnose

    N_orig = length(ir_original);
    t_orig = (0:N_orig-1) / fs * 1000;  % Zeit in ms

    % Zeitlimit: Ende der Truncation + 10% Puffer, maximal 500ms
    t_end = t_orig(end_idx);
    t_max = min(500, t_end * 1.1);  % Ende + 10% Puffer, maximal 500ms

    fig = figure('Visible','off','Position', [100, 100, 1400, 800]);

    % Subplot 1: Originale IR (Amplitude)
    subplot(3,1,1);
    plot(t_orig, ir_original, 'b-', 'LineWidth', 1);
    hold on;
    % Markiere Start und Ende
    xline(t_orig(start_idx), 'g--', 'LineWidth', 2, 'Label', 'Start');
    xline(t_orig(end_idx), 'r--', 'LineWidth', 2, 'Label', 'Ende');
    hold off;
    grid on;
    xlabel('Zeit [ms]');
    ylabel('Amplitude');
    title(sprintf('%s - Position %02d: Energie %.2f%%, SNR %.1f dB, DR %.1f dB', ...
                  variantName, pos, E_ratio, SNR_dB, dynamic_range_dB));
    xlim([0 t_max]);

    % Subplot 2: Originale IR (dB)
    subplot(3,1,2);
    ir_abs = abs(ir_original);
    ir_dB = 20*log10(ir_abs / max(ir_abs) + eps);
    plot(t_orig, ir_dB, 'b-', 'LineWidth', 1);
    hold on;
    xline(t_orig(start_idx), 'g--', 'LineWidth', 2, 'Label', 'Start');
    xline(t_orig(end_idx), 'r--', 'LineWidth', 2, 'Label', 'Ende');
    yline(-60, 'k:', 'LineWidth', 1.5, 'Label', '-60 dB');
    hold off;
    grid on;
    xlabel('Zeit [ms]');
    ylabel('Pegel [dB]');
    title('Pegel relativ zum Maximum');
    xlim([0 t_max]);
    ylim([-100 5]);

    % Subplot 3: Truncierte IR (Amplitude)
    subplot(3,1,3);
    N_trunc = length(ir_trunc);
    t_trunc = (0:N_trunc-1) / fs * 1000;  % Zeit in ms
    plot(t_trunc, ir_trunc, 'r-', 'LineWidth', 1);
    grid on;
    xlabel('Zeit [ms]');
    ylabel('Amplitude');
    title(sprintf('Truncierte Impulsantwort (Länge: %d → %d Samples)', N_orig, N_trunc));
    xlim([0 max(t_trunc)]);

    % Speichern
    filename_png = fullfile('Plots', sprintf('%s_Pos_%02d_E%.1f.png', variantName, pos, E_ratio));
    filename_fig = fullfile('Plots', sprintf('%s_Pos_%02d_E%.1f.fig', variantName, pos, E_ratio));

    saveas(fig, filename_png);
    savefig(fig, filename_fig);
    close(fig);

    fprintf('  → Diagnose-Plot erstellt: %s\n', filename_png);
end

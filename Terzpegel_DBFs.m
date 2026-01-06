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
selectedVariants = { 'Variante_2'};  % z.B. {'Variante_1', 'Variante_2'} oder {} für Auto
excludePattern = '_alt';  % Ausschlussmuster für Varianten

% Positionen (Messpunkte) auswählen
selectedPositions = 1:14;   % z.B. [1 3 5] oder 1:14

% Ausgabeordner
outputPlotDir = 'Plots';
outputExcelDir = 'Excel';

% Plot-Modus: 'absolute' (Pegel in dBFS relativ zur globalen Referenz)
%            'difference' (je Messung: Pegel minus Mean-Pegel über Positionen für diese Variante)
plotMode = 'absolute';   % oder 'difference'

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

        fprintf('  Position %02d ausgewertet\n', pos);

        N = length(ir);
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

    % Optional: Differenzmodus
    L_plot = L_dBFS;
    if strcmpi(plotMode,'difference')
        L_plot = L_plot - L_mean_dBFS;
    end

    %% ---------------- Excel-Export ----------------
    rowNames = compose('Pos_%02d', selectedPositions);
    colNames = compose('F%.0f', f_terz);

    T_LdBFS = array2table(L_plot, 'RowNames', rowNames, 'VariableNames', colNames);

    excelFile = fullfile(outputExcelDir, sprintf('%s_Terzpegel_dBFS.xlsx', variantName));
    try
        writetable(T_LdBFS, excelFile, 'WriteRowNames', true);
        fprintf('Excel geschrieben: %s\n', excelFile);
    catch ME
        warning(ME.identifier, 'Konnte Excel nicht schreiben: %s', ME.message);
    end

    %% ---------------- Plots ----------------
    % Y-Achsen-Grenzen (global für diese Variante)
    validData = L_plot(~isnan(L_plot));
    if ~isempty(validData)
        y_range_terz = [floor(min(validData)/10)*10, ceil(max(validData)/10)*10];
    else
        y_range_terz = [-100 0];
    end

    % Plot-Einstellungen (konstant für alle Plots)
    xtick_vals = [500 1000 2000 5000 10000 20000 50000 100000];
    xtick_labels = {'500','1k','2k','5k','10k','20k','50k','100k'};
    f_lim = [min(f_terz) max(f_terz)];

    % Dateinamen-Suffix basierend auf Modus
    if strcmpi(plotMode, 'difference')
        modeSuffix = '_diff';
    else
        modeSuffix = '';
    end

    % Terzpegel-Plots
    for pi = 1:nPos
        pos = selectedPositions(pi);
        if all(isnan(L_plot(pi,:))), continue; end

        fig = figure('Visible','off','Position',[100,100,1000,500]);
        stairs(f_terz, L_plot(pi,:), 'LineWidth',2, 'Color',[0 0.4470 0.7410]);
        grid on; set(gca,'XScale','log');
        xlabel('Frequenz [Hz]'); ylabel('Pegel [dB]');
        title(sprintf('%s - Position %02d', variantName, pos));
        xlim(f_lim);
        ylim(y_range_terz);
        xticks(xtick_vals);
        xticklabels(xtick_labels);

        filename = fullfile(outputPlotDir, sprintf('Terzpegel_%s_Pos_%02d%s.png', variantName, pos, modeSuffix));
        saveas(fig, filename);
        saveas(fig, strrep(filename,'.png','.fig'));
        close(fig);
    end

    % Übertragungsfunktionen
    validH = cellfun(@(x) ~isempty(x), H_all);
    if any(validH)
        allH = cell2mat(H_all(validH)');
        y_range_H = [floor(min(allH(:))/10)*10, ceil(max(allH(:))/10)*10];
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

    fprintf('%s abgeschlossen: Excel + Plots gespeichert.\n', variantName);
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

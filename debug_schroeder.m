%% ============================================================
%  Debug: Schröder-Kurven visualisieren
%  Zeigt Schröder-Integration und T20-Regressionsbereiche
% ============================================================

clear;
clc;

%% ---------------- Einstellungen ----------------
fs = 500e3;        % Abtastrate [Hz]

% Normgerechte 1/3-Oktav-Mittenfrequenzen (IEC 61260)
f_terz = [ ...
    560 625 710 800 900 1000 1120 1250 1400 1600 1800 2000 ...
    2240 2500 2800 3150 3550 4000 4500 5000 5600 6300 7100 8000 ...
    9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 ...
    28000 31500 35500 40000 45000 50000 56000 63000 71000 80000 ...
    90000 100000 112000 126000 ];

filterOrder = 4;

% --- ANPASSEN: Welche Variante und Position soll geplottet werden? ---
variantName = 'Variante_4';
messPos = 1;  % Position 1-14

% --- ANPASSEN: Welche Terzbänder sollen geplottet werden? ---
% z.B. [10 20 30] für 3 ausgewählte Bänder, oder 1:length(f_terz) für alle
plotBands = [10 15 20 25 30 35];  % 6 Beispiel-Bänder

%% ---------------- Datei laden ----------------
filename = sprintf('%s/%s,Pos_%d.mat', variantName, variantName, messPos);

if ~isfile(filename)
    error('Datei nicht gefunden: %s', filename);
end

data = load(filename, 'RIR');

if ~isfield(data, 'RIR')
    error('Datei %s enthält keinen Vektor "RIR".', filename);
end

ir = data.RIR(:);
t = (0:length(ir)-1)' / fs;

fprintf('Datei geladen: %s\n', filename);
fprintf('Dauer: %.3f s, Samples: %d\n\n', t(end), length(ir));

%% ---------------- Plot-Vorbereitung ----------------
nPlots = length(plotBands);
nCols = 2;
nRows = ceil(nPlots / nCols);

figure('Position', [100 100 1200 800]);

%% ---------------- Schröder-Kurven plotten ----------------
for idx = 1:nPlots

    k = plotBands(idx);

    if k > length(f_terz)
        warning('Terzband-Index %d außerhalb des Bereichs.', k);
        continue;
    end

    % --- Terzband-Grenzen ---
    f1 = f_terz(k) / 2^(1/6);
    f2 = f_terz(k) * 2^(1/6);
    Wn = [f1 f2] / (fs/2);

    if Wn(2) >= 1
        fprintf('Terzband %d Hz: Oberhalb Nyquist, übersprungen.\n', f_terz(k));
        continue;
    end

    % --- Terzbandfilter ---
    [b, a] = butter(filterOrder/2, Wn, 'bandpass');
    ir_filt = filtfilt(b, a, ir);

    % --- Schroeder-Integration ---
    E = flipud(cumsum(flipud(ir_filt.^2)));
    E = E / max(E);
    E_dB = 10*log10(E + eps);

    % --- T20-Bereich ---
    idx_t20 = find(E_dB <= -5 & E_dB >= -25);

    % --- Lineare Regression ---
    RT60 = NaN;
    if numel(idx_t20) >= 10
        p = polyfit(t(idx_t20), E_dB(idx_t20), 1);
        RT60 = -60 / p(1);

        % Regressionsgerade
        t_fit = [t(idx_t20(1)); t(idx_t20(end))];
        E_fit = polyval(p, t_fit);
    end

    % --- Subplot ---
    subplot(nRows, nCols, idx);

    % Schröder-Kurve
    plot(t, E_dB, 'b-', 'LineWidth', 1.5); hold on;

    % T20-Bereich markieren
    if numel(idx_t20) >= 10
        plot(t(idx_t20), E_dB(idx_t20), 'r.', 'MarkerSize', 8);
        plot(t_fit, E_fit, 'g--', 'LineWidth', 2);
    end

    % Horizontale Linien bei -5, -25, -60 dB
    yline(-5, 'k:', 'Alpha', 0.3);
    yline(-25, 'k:', 'Alpha', 0.3);
    yline(-60, 'k:', 'Alpha', 0.3);

    % Beschriftung
    ylim([-70 5]);
    xlim([0 t(end)]);
    grid on;
    xlabel('Zeit [s]');
    ylabel('Energie [dB]');
    title(sprintf('f = %d Hz | RT60 = %.2f s', f_terz(k), RT60));
    legend('Schröder-Kurve', 'T20-Bereich', 'Regression', ...
           'Location', 'northeast');

    fprintf('Terzband %5d Hz: RT60 = %.3f s, T20-Punkte = %d\n', ...
            f_terz(k), RT60, numel(idx_t20));
end

sgtitle(sprintf('%s, Position %d - Schröder-Kurven Debug', ...
                variantName, messPos), 'FontSize', 14, 'FontWeight', 'bold');

fprintf('\n✓ Plots erstellt!\n');

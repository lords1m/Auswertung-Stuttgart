%% ============================================================
%  Optimierte Energy Decay Curve (EDC) Berechnung
%  Methode: Schroeder-Integration mit Rauschkompensation
%% ============================================================

clear;
clc;
close all;

%% ---------------- 1. Einstellungen ----------------
fs = 500e3;        % Abtastrate [Hz] (Bitte prüfen, ob 500kHz korrekt ist!)
filename = 'Variante_4/Variante_4,Pos_1.mat';

%% ---------------- 2. Datei laden ----------------
fprintf('Lade Datei: %s\n', filename);
try
    data = load(filename, 'RIR');
    ir_raw = data.RIR(:);
catch
    error('Datei konnte nicht geladen werden oder enthält kein Feld "RIR".');
end

%% ---------------- 3. Vorverarbeitung (Preprocessing) ----------------
% 3.1 Startpunkt finden (Direktschall-Peak)
[~, maxIdx] = max(abs(ir_raw));
% Wir schneiden alles vor dem Peak ab, um die Totzeit zu eliminieren
ir = ir_raw(maxIdx:end);

% 3.2 Zeitvektoren erstellen
t_raw = (0:length(ir_raw)-1)' / fs;
t_edc = (0:length(ir)-1)' / fs;

%% ---------------- 4. EDC Berechnung (Schroeder) ----------------
% Quadrat der Impulsantwort (Energie)
ir_sq = ir.^2;

% --- RAUSCHKOMPENSATION ---
% Wir schätzen das Rauschen aus den letzten 10% der Aufnahme
noise_floor_est = mean(ir_sq(round(end*0.9):end));

% Subtraktionsmethode: Rauschen abziehen, um die EDC zu begradigen
ir_sq_compensated = ir_sq - noise_floor_est;
ir_sq_compensated(ir_sq_compensated < 0) = 0; % Keine negativen Energiewerte

% Rückwärts-Integration (Schroeder)
E = flipud(cumsum(flipud(ir_sq_compensated)));

% Normalisierung (Max = 0 dB)
E_norm = E / max(E);
E_dB = 10 * log10(E_norm + eps);

fprintf('EDC berechnet und Rauschen kompensiert.\n');

%% ---------------- 6. RT60 Berechnung (T20/T10 Methode) ----------------

% Ziel-Pegel für T20 (von -5 bis -25 dB)
% Da dein Rauschen bei -30dB liegt, ist T20 riskant. T10 ist sicherer.
target_start = -5;
target_end   = -15; % Wir nehmen T10 (von -5 bis -15), weil das SNR klein ist

% Indizes finden, wo die EDC diese Pegel kreuzt
idx_start = find(E_dB <= target_start, 1, 'first');
idx_end   = find(E_dB <= target_end, 1, 'first');

if isempty(idx_start) || isempty(idx_end)
    warning('Signal fällt nicht tief genug für die Berechnung!');
    rt60_val = NaN;
else
    % Zeitpunkte extrahieren
    t_segment = t_edc(idx_start:idx_end);
    E_segment = E_dB(idx_start:idx_end);
    
    % Lineare Regression: E_dB = m * t + b
    p = polyfit(t_segment, E_segment, 1);
    slope = p(1); % Abfall pro Sekunde in dB
    
    % RT60 ist die Zeit, in der das Signal um 60 dB fallen würde
    % Formel: -60 / Steigung
    rt60_val = -60 / slope;
    
    fprintf('Berechnetes RT60 (basiert auf T10): %.3f Sekunden\n', rt60_val);
    
    % Regressionsgerade für den Plot berechnen
    t_fit = linspace(0, rt60_val, 100);
    E_fit = polyval(p, t_fit);
end

%% ---------------- 7. Plot mit RT60-Linie ----------------
figure('Position', [100 100 800 500]);
plot(t_edc, E_dB, 'LineWidth', 1.5, 'DisplayName', 'EDC Kurve');
hold on;

if ~isnan(rt60_val)
    plot(t_fit, E_fit, 'k--', 'LineWidth', 2, 'DisplayName', sprintf('Fit (RT60 = %.2fs)', rt60_val));
    % Markierung der Messpunkte
    plot(t_edc(idx_start), E_dB(idx_start), 'go', 'MarkerFaceColor', 'g');
    plot(t_edc(idx_end), E_dB(idx_end), 'ro', 'MarkerFaceColor', 'r');
end

grid on;
ylim([-60 5]);
xlim([0 max(t_edc(1:round(end/2)))]); % Zoom auf den relevanten Teil
xlabel('Zeit [s]');
ylabel('Energie [dB]');
title(['RT60 Schätzung: ' num2str(rt60_val, '%.3f') ' s']);
legend show;

%% ---------------- 5. Visualisierung ----------------
figure('Position', [100 100 1000 700], 'Name', 'Akustische Analyse');

% Subplot 1: Rohdaten (Logarithmisch zur Rauschanalyse)
subplot(2,1,1);
ir_log = 20 * log10(abs(ir_raw) / max(abs(ir_raw)) + eps);
plot(t_raw, ir_log, 'Color', [0.5 0.5 0.5]); % Grau für Rohdaten
hold on;
% Markierung des Startpunkts
plot(t_raw(maxIdx), 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2); 
grid on;
title('Logarithmische Impulsantwort (Rohdaten)');
ylabel('Pegel [dB]');
xlabel('Zeit [s]');
legend('RIR', 'Detektierter Start (Peak)');
ylim([-100 5]);

% Subplot 2: Energy Decay Curve
subplot(2,1,2);
plot(t_edc, E_dB, 'r-', 'LineWidth', 2);
grid on;
title('Energy Decay Curve (Schroeder-Integration)');
ylabel('Energie [dB]');
xlabel('Zeit [s]');
ylim([-60 5]); % Fokus auf die relevanten 60dB Abfall
xlim([0 max(t_edc)]);

% Titel über das ganze Fenster
sgtitle(['Analyse: ' strrep(filename, '_', '\_')]);

fprintf('Fertig. Wenn die Kurve immer noch flach ist, ist das SNR zu gering.\n');
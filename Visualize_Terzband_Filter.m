%% ============================================================
%  Visualisierung des Terzband-Filters
%  Zeigt die Filterkurve und gefilterte Impulsantwort
% ============================================================

clear;
clc;

%% ---------------- Einstellungen (HIER ANPASSEN) ----------------

% Welche Variante und Position?
selectedVariant = 'Variante_3';
selectedPosition = 9;

% Welche Terzband-Mittenfrequenz visualisieren? [Hz]
selectedFrequency = 8000;  % z.B. 1000, 10000, 50000

% Datenordner
dataDir = 'data';

% Abtastrate
fs = 500e3;  % 500 kHz

% Filterordnung (wie im Hauptskript)
filterOrder = 4;

%% ---------------- Lade Impulsantwort ----------------

% Dateiname konstruieren
filename = fullfile(dataDir, sprintf('%s,Pos_%d.mat', selectedVariant, selectedPosition));

if ~exist(filename, 'file')
    error('Datei nicht gefunden: %s', filename);
end

fprintf('Lade Datei: %s\n', filename);
S = load(filename);

% Extrahiere Impulsantwort (wie im Hauptskript)
if isfield(S, 'RIR')
    ir = S.RIR(:);
elseif isfield(S, 'IR')
    ir = S.IR(:);
elseif isfield(S, 'ir')
    ir = S.ir(:);
else
    error('Keine Impulsantwort gefunden. Erwartete Felder: RIR, IR, oder ir');
end

N_original = length(ir);
fprintf('Impulsantwort geladen: %d Samples (%.2f ms)\n', N_original, N_original/fs*1000);

%% ---------------- Terzband-Parameter ----------------

% Normgerechte 1/3-Oktav-Mittenfrequenzen (IEC 61260)
f_terz = double([ ...
    560 625 710 800 900 1000 1120 1250 1400 1600 1800 2000 ...
    2240 2500 2800 3150 3550 4000 4500 5000 5600 6300 7100 8000 ...
    9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 ...
    28000 31500 35500 40000 45000 50000 56000 63000 71000 80000 ...
    90000 100000 112000 126000 ]);

% Finde nächste Terzband-Mittenfrequenz
[~, freq_idx] = min(abs(f_terz - selectedFrequency));
f_center = f_terz(freq_idx);

fprintf('\nGewählte Mittenfrequenz: %.0f Hz\n', f_center);

% Terzband-Grenzen berechnen
f_lower = f_center / 2^(1/6);
f_upper = f_center * 2^(1/6);

fprintf('Terzband-Grenzen: %.1f Hz - %.1f Hz\n', f_lower, f_upper);

% Normierte Grenzfrequenzen
Wn = [f_lower f_upper] / (fs/2);

if Wn(2) >= 1
    error('Obere Grenzfrequenz (%.0f Hz) liegt über Nyquist-Frequenz (%.0f Hz)!', f_upper, fs/2);
end

%% ---------------- Butterworth-Bandpass-Filter ----------------

% Filter-Design (wie im Hauptskript)
[b, a] = butter(filterOrder/2, Wn, 'bandpass');

fprintf('Filter erstellt: Butterworth Bandpass, Ordnung %d\n', filterOrder);

%% ---------------- Filter-Frequenzgang ----------------

% Berechne Frequenzgang des Filters
[H_filter, f_filter] = freqz(b, a, 4096, fs);
H_filter_dB = 20*log10(abs(H_filter));

%% ---------------- Filtere Impulsantwort ----------------

% Filtere mit filtfilt (zero-phase)
ir_filtered = filtfilt(b, a, ir);

% Zeitvektor
t = (0:N_original-1) / fs * 1000;  % Zeit in ms

%% ---------------- FFT der Original- und gefilterten IR ----------------

% FFT
N_fft = 2^nextpow2(N_original * 2);  % Zero-padding für bessere Frequenzauflösung
IR_fft = fft(ir, N_fft);
IR_filtered_fft = fft(ir_filtered, N_fft);

% Frequenzvektor
freq = (0:N_fft-1) * (fs / N_fft);
nHalf = floor(N_fft/2) + 1;
freq = freq(1:nHalf);

% Magnitude in dB
H_original = 20*log10(abs(IR_fft(1:nHalf)) + eps);
H_filtered = 20*log10(abs(IR_filtered_fft(1:nHalf)) + eps);

% Normiere auf Maximum
H_original = H_original - max(H_original);
H_filtered = H_filtered - max(H_filtered);

%% ---------------- Plots ----------------

figure('Position', [100, 100, 1400, 900]);

% Subplot 1: Filter-Frequenzgang
subplot(3,2,1);
semilogx(f_filter, H_filter_dB, 'b-', 'LineWidth', 2);
hold on;
xline(f_lower, 'r--', 'LineWidth', 1.5, 'Label', sprintf('f_l=%.1f Hz', f_lower));
xline(f_upper, 'r--', 'LineWidth', 1.5, 'Label', sprintf('f_u=%.1f Hz', f_upper));
xline(f_center, 'g-', 'LineWidth', 2, 'Label', sprintf('f_c=%.0f Hz', f_center));
yline(-3, 'k:', 'LineWidth', 1, 'Label', '-3 dB');
hold off;
grid on;
xlabel('Frequenz [Hz]');
ylabel('Dämpfung [dB]');
title(sprintf('Butterworth Bandpass Filter (Ordnung %d)', filterOrder));
xlim([f_lower/2, f_upper*2]);
ylim([-60 5]);

% Subplot 2: Zoom auf Durchlassbereich
subplot(3,2,2);
semilogx(f_filter, H_filter_dB, 'b-', 'LineWidth', 2);
hold on;
xline(f_lower, 'r--', 'LineWidth', 1.5);
xline(f_upper, 'r--', 'LineWidth', 1.5);
xline(f_center, 'g-', 'LineWidth', 2);
yline(-3, 'k:', 'LineWidth', 1, 'Label', '-3 dB');
hold off;
grid on;
xlabel('Frequenz [Hz]');
ylabel('Dämpfung [dB]');
title('Filter (Zoom auf Durchlassbereich)');
xlim([f_lower/1.5, f_upper*1.5]);
ylim([-6 1]);

% Subplot 3: Original Impulsantwort (Zeitbereich)
subplot(3,2,3);
plot(t, ir, 'b-', 'LineWidth', 1);
grid on;
xlabel('Zeit [ms]');
ylabel('Amplitude');
title(sprintf('Original Impulsantwort - %s, Position %d', selectedVariant, selectedPosition));
xlim([0 min(100, max(t))]);  % Erste 100ms oder weniger

% Subplot 4: Gefilterte Impulsantwort (Zeitbereich)
subplot(3,2,4);
plot(t, ir_filtered, 'r-', 'LineWidth', 1);
grid on;
xlabel('Zeit [ms]');
ylabel('Amplitude');
title(sprintf('Gefilterte Impulsantwort (f_c = %.0f Hz)', f_center));
xlim([0 min(100, max(t))]);

% Subplot 5: Original Spektrum
subplot(3,2,5);
semilogx(freq, H_original, 'b-', 'LineWidth', 1);
hold on;
xline(f_lower, 'r--', 'LineWidth', 1.5);
xline(f_upper, 'r--', 'LineWidth', 1.5);
xline(f_center, 'g-', 'LineWidth', 2);
hold off;
grid on;
xlabel('Frequenz [Hz]');
ylabel('Pegel [dB rel. Max]');
title('Original Spektrum (FFT)');
xlim([4000 60000]);
ylim([-80 5]);

% Subplot 6: Gefiltertes Spektrum
subplot(3,2,6);
semilogx(freq, H_filtered, 'r-', 'LineWidth', 1.5);
hold on;
xline(f_lower, 'r--', 'LineWidth', 1.5, 'Label', sprintf('%.1f Hz', f_lower));
xline(f_upper, 'r--', 'LineWidth', 1.5, 'Label', sprintf('%.1f Hz', f_upper));
xline(f_center, 'g-', 'LineWidth', 2, 'Label', sprintf('%.0f Hz', f_center));
hold off;
grid on;
xlabel('Frequenz [Hz]');
ylabel('Pegel [dB rel. Max]');
title('Gefiltertes Spektrum (FFT)');
xlim([4000 60000]);
ylim([-80 5]);

% Gesamttitel
sgtitle(sprintf('Terzband-Filter Visualisierung: %s, Pos %d, f_c = %.0f Hz', ...
                selectedVariant, selectedPosition, f_center), ...
        'FontSize', 14, 'FontWeight', 'bold');

%% ---------------- Speichern ----------------

outputDir = 'Plots';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

filename_out = fullfile(outputDir, sprintf('Filter_%s_Pos%02d_f%.0f.png', ...
                                            selectedVariant, selectedPosition, f_center));
saveas(gcf, filename_out);
saveas(gcf, strrep(filename_out, '.png', '.fig'));

fprintf('\n========================================\n');
fprintf('Filter-Visualisierung gespeichert:\n');
fprintf('  %s\n', filename_out);
fprintf('========================================\n');

%% ---------------- Zusätzliche Informationen ----------------

fprintf('\nFilter-Informationen:\n');
fprintf('  Mittenfrequenz:     %.0f Hz\n', f_center);
fprintf('  Untere Grenze:      %.1f Hz\n', f_lower);
fprintf('  Obere Grenze:       %.1f Hz\n', f_upper);
fprintf('  Bandbreite:         %.1f Hz\n', f_upper - f_lower);
fprintf('  Relative Bandbreite: %.2f%%\n', (f_upper - f_lower) / f_center * 100);
fprintf('  Filterordnung:      %d\n', filterOrder);
fprintf('  Abtastrate:         %.0f kHz\n', fs/1000);
fprintf('\n');

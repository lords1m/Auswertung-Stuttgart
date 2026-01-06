%% Test der EDC-Berechnung (Schröder-Integration)
clear; clc;

% Lade eine Beispiel-RIR
fs = 500e3;
data = load('Variante_4/Variante_4,Pos_1.mat', 'RIR');
ir = data.RIR(:);
t = (0:length(ir)-1)' / fs;

% Test-Frequenz: 1000 Hz
f_center = 1000;
f1 = f_center / 2^(1/6);
f2 = f_center * 2^(1/6);
Wn = [f1 f2] / (fs/2);

% Terzbandfilter
[b, a] = butter(2, Wn, 'bandpass');
ir_filt = filtfilt(b, a, ir);

%% Verschiedene EDC-Berechnungen zum Vergleich

% Methode 1: Original (flip-cumsum-flip)
ir_sq = ir_filt.^2;
E1 = flip(cumsum(flip(ir_sq)));
E1 = E1 / max(E1);
E1_dB = 10*log10(E1 + eps);

% Methode 2: flipud-cumsum-flipud
E2 = flipud(cumsum(flipud(ir_sq)));
E2 = E2 / max(E2);
E2_dB = 10*log10(E2 + eps);

% Methode 3: Rückwärts-Integration (korrekt nach Schroeder 1965)
E3 = zeros(size(ir_sq));
for i = 1:length(ir_sq)
    E3(i) = sum(ir_sq(i:end));
end
E3 = E3 / max(E3);
E3_dB = 10*log10(E3 + eps);

% Methode 4: cumsum von hinten nach vorne
E4 = cumsum(ir_sq(end:-1:1));
E4 = E4(end:-1:1);
E4 = E4 / max(E4);
E4_dB = 10*log10(E4 + eps);

%% Plot zum Vergleich
figure('Position', [100 100 1400 800]);

subplot(2,2,1);
plot(t*1000, E1_dB, 'b-', 'LineWidth', 1.5);
grid on; xlabel('Zeit [ms]'); ylabel('Energie [dB]');
title('Methode 1: flip(cumsum(flip(ir^2)))');
ylim([-60 5]); xlim([0 1000]);

subplot(2,2,2);
plot(t*1000, E2_dB, 'r-', 'LineWidth', 1.5);
grid on; xlabel('Zeit [ms]'); ylabel('Energie [dB]');
title('Methode 2: flipud(cumsum(flipud(ir^2)))');
ylim([-60 5]); xlim([0 1000]);

subplot(2,2,3);
plot(t*1000, E3_dB, 'g-', 'LineWidth', 1.5);
grid on; xlabel('Zeit [ms]'); ylabel('Energie [dB]');
title('Methode 3: Explizite Rückwärts-Integration (Referenz)');
ylim([-60 5]); xlim([0 1000]);

subplot(2,2,4);
plot(t*1000, E4_dB, 'm-', 'LineWidth', 1.5);
grid on; xlabel('Zeit [ms]'); ylabel('Energie [dB]');
title('Methode 4: cumsum(ir(end:-1:1)^2) rückgekehrt');
ylim([-60 5]); xlim([0 1000]);

sgtitle(sprintf('EDC-Vergleich: 1000 Hz Terzband'));

% Vergleich der Ergebnisse
fprintf('Vergleich der Methoden:\n');
fprintf('Methode 1 == Methode 2: %d\n', all(abs(E1_dB - E2_dB) < 1e-10));
fprintf('Methode 1 == Methode 3: %d\n', all(abs(E1_dB - E3_dB) < 1e-10));
fprintf('Methode 3 == Methode 4: %d\n', all(abs(E3_dB - E4_dB) < 1e-10));

% Speichere Plot
saveas(gcf, 'EDC_Vergleich.png');
fprintf('\nPlot gespeichert als EDC_Vergleich.png\n');

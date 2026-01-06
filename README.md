# Auswertung Stuttgart - MATLAB Skripte für Raumimpulsantwort-Analyse

Dieses Repository enthält MATLAB-Skripte zur Analyse von Raumimpulsantworten (RIR) mit Fokus auf Terzpegel-Auswertung und Nachhallzeit-Berechnung.

Variante 1 = Ohne Elemente, Ohne Deckel
Variante 2 = Ohne Elemente, Mit Deckel
Variante 3 = Mit Absorber, Ohne Deckel
Variante 4 = Mit Balkonen, Ohne Deckel
---

Alle Plots stellen den bereich 4k - 60k dar. 

## Erforderliche Ordnerstruktur

```
Auswertung-Stuttgart/
├── data/                          # Eingabedaten (muss vorhanden sein!)
│   ├── Variante_1_neu,Pos_1.mat
│   ├── Variante_1_neu,Pos_2.mat
│   ├── ...
│   ├── Variante_2,Pos_1.mat
│   └── ...
├── Plots/                         # Ausgabe: Alle Plots (wird automatisch erstellt)
├── Excel/                         # Ausgabe: Excel-Dateien (wird automatisch erstellt)
├── Terzpegel_DBFs.m              # Hauptskript: Terzpegel-Analyse
├── RT_60.m                        # Nachhallzeit-Berechnung
├── Visualize_Terzband_Filter.m   # Filter-Visualisierung
├── Terzpegel_alle_Varianten.m    # (Optional) Varianten-Vergleich
└── README.md                      # Diese Datei
```

## Skript-Übersicht

### 1. Terzpegel_DBFs.m - Terzpegel-Analyse (Hauptskript)

**Zweck:**
Berechnet Terzpegel (1/3-Oktavband-Analyse) für Raumimpulsantworten aller Varianten und Positionen.

**Funktionen:**
- Lundeby-Truncation zur Rauschunterdrückung
- FFT-basierte Spektralanalyse
- Terzpegel-Berechnung (IEC 61260 normgerecht)
- dBFS-Normierung relativ zum globalen Maximum
- Mittelwert- und Summen-Terzpegel
- Mittelwert-Spektrum über alle Positionen
- SNR und Dynamic Range Logging
- Automatische Diagnose-Plots bei niedriger Energie (<60%)

**Einstellungen (Zeilen 22-31):**
```matlab
dataDir = 'data';                                    % Verzeichnis des Datenordner
selectedVariants = {'Variante_1_neu', 'Variante_2'}; % Diese Variantten sollen ausgewertet werden
selectedPositions = 1:14;                            % Diese Positionen sollen ausgewertet werden
outputPlotDir = 'Plots';                             % Hier sollen die Plots rein
outputExcelDir = 'Excel';                            % Hier sollen die Excel-Arbeitsmappen rein
plotMode = 'absolute';                               % Die Werte sollen absolut oder im Verhältnis zum maximalwert berechnet werden
```

**Ausgabe:**
- **Excel:** `Excel/[Variante]_Terzpegel_dBFS.xlsx` (14 Positionen × 48 Terzbänder)
- **Plots:**
  - `Plots/Terzpegel_[Variante]_Pos_[XX].png` (14 Plots pro Variante)
  - `Plots/Uebertragungsfunktion_[Variante]_Pos_[XX].png` (FFT-Spektrum)
  - `Plots/Terzpegel_[Variante]_Summe.png` (Energetische Summe)
  - `Plots/Spektrum_[Variante]_Mittelwert.png` (Gemitteltes Spektrum)
  - `Plots/[Variante]_Pos_[XX]_E[XX.X].png` (Diagnose bei <60% Energie)


2. **RT_60.m** - Nachhallzeit-Berechnung

Berechnet die Nachhallzeit RT60 für alle Varianten und Positionen mittels Schroeder-Integration und T20-Methode.

**Funktionen:**
- Lundeby-Truncation
- Terzband-Filterung (IEC 61260)
- Schroeder-Integration für Energieabfall
- T20-Bereich (-5 dB bis -25 dB) mit linearer Regression
- RT60-Extrapolation (-60 dB)
- Mittelwert über ausgewählte Positionen (5,6,7,9,10,11,13,14)
- Varianten-Vergleich

**Einstellungen (Zeilen 11-30):**
```matlab
nMess = 14;                    % Anzahl Messpositionen
fs = 500e3;                    % Abtastrate [Hz]
startFolder = 'Variante_1_neu'; % Ab welcher Variante starten?
positions_to_average = [5, 6, 7, 9, 10, 11, 13, 14]; % Positionen für Mittelwertbildung
```

### 3. **Visualize_Terzband_Filter.m** - Filter-Visualisierung

**Zweck:**
Visualisiert den Butterworth-Bandpass-Filter für eine spezifische Variante, Position und Terzband-Mittenfrequenz.

**Funktionen:**
- Filter-Frequenzgang (breit + Zoom)
- Original vs. gefilterte Impulsantwort (Zeit + Frequenz)
- Markierungen für Mittenfrequenz, Bandgrenzen und -3 dB Punkt
- 6 Subplots für vollständige Analyse

**Einstellungen (Zeilen 11-16):**
```matlab
selectedVariant = 'Variante_2';   % Welche Variante?
selectedPosition = 1;              % Welche Position?
selectedFrequency = 10000;         % Welche Frequenz? [Hz]
```

**Ausgabe:**
- **Plot:** `Plots/Filter_[Variante]_Pos[XX]_f[XXXX].png/.fig`

### 4. **Terzpegel_alle_Varianten.m** (Optional)

**Zweck:**
Vergleicht Terzpegel mehrerer Varianten in einem Plot.

**Status:** Vorhanden, aber möglicherweise veraltet. Verwende stattdessen die Vergleichsfunktion in `Terzpegel_DBFs.m`.

## Dateiformat-Anforderungen

### MAT-Dateien

Jede `.mat`-Datei muss eine Raumimpulsantwort enthalten:

```matlab
% Akzeptierte Variablennamen (in Priorität):
RIR   % Bevorzugt
IR    % Alternative
ir    % Alternative
```

**Dateinamen-Format:**
```
[VariantenName],Pos_[Nummer].mat
```

**Beispiele:**
- `Variante_1_neu,Pos_1.mat`
- `Variante_2,Pos_15.mat`
- `Variante_3_alt,Pos_7.mat` (wird mit `excludePattern` gefiltert)

---

## Terzband-Frequenzen (IEC 61260)

Die Skripte verwenden 48 normgerechte 1/3-Oktavband-Mittenfrequenzen:

```
560 Hz, 625 Hz, 710 Hz, ..., 100 kHz, 112 kHz, 126 kHz
```

**Bandgrenzen:**
- Untere Grenze: `f_center / 2^(1/6)`
- Obere Grenze: `f_center * 2^(1/6)`

---

## Wichtige Berechnungen

### Lundeby-Truncation
Entfernt Rauschen am Anfang und Ende der Impulsantwort:
- **Startpunkt:** Erster Peak über 5% des Maximums
- **Endpunkt:** Letzter Peak über 10× Rauschpegel + 20% Sicherheitsmarge
- **Rauschpegel:** Mittelwert + 3σ der letzten 10% der Samples

### SNR (Signal-to-Noise Ratio)
```matlab
SNR [dB] = 20 × log10(signal_rms / noise_rms)
```
- `signal_rms`: RMS der truncierten Impulsantwort
- `noise_rms`: Standardabweichung der letzten 10%

### Dynamic Range
```matlab
DR [dB] = 20 × log10(max_amplitude / noise_level)
```

### Terzpegel-Mittelwert (energetisch)
```matlab
L_mean [dB] = 10 × log10( Σ(10^(L_i/10)) / N )
```

### Summen-Terzpegel (energetische Addition)
```matlab
L_sum [dB] = 10 × log10( Σ(10^(L_i/10)) )
```

---

## Systemanforderungen

- **MATLAB:** R2019b oder neuer (wegen `compose()`, `xline()`)
- **Toolboxen:**
  - Signal Processing Toolbox (für `butter()`, `filtfilt()`, `freqz()`)
  - (Optional) Statistics Toolbox (für robustere Auswertungen)

---

## Schnellstart

1. **Ordner erstellen:**
   ```bash
   mkdir data Plots Excel
   ```

2. **Daten kopieren:**
   Alle `.mat`-Dateien in `data/` Ordner kopieren.

3. **Skript anpassen:**
   ```matlab
   % In Terzpegel_DBFs.m:
   selectedVariants = {'Variante_1_neu', 'Variante_2'};
   selectedPositions = 1:14;
   ```

4. **Ausführen:**
   ```matlab
   run('Terzpegel_DBFs.m')
   ```

5. **Ergebnisse prüfen:**
   - Excel-Dateien in `Excel/`
   - Plots in `Plots/`

---

## Hinweise

### Energie-Warnung (<60%)
Wenn eine Impulsantwort weniger als 60% Energie nach Truncation behält, wird automatisch ein Diagnose-Plot erstellt:
```
Plots/[Variante]_Pos_[XX]_E[XX.X].png
```

**Grenzwert anpassen:**
- Datei: `Terzpegel_DBFs.m`
- Zeile: 198 - `if E_ratio < 60` (60% ändern für andere Schwelle)

### Plots werden nicht angezeigt
Alle Plots verwenden `'Visible','off'` und werden direkt gespeichert, ohne Fenster zu öffnen.

**Plots sichtbar machen:**
- Datei: `Terzpegel_DBFs.m`
- Zeilen: 277, 314, 334, 371, 477 - `'Visible','off'` auf `'Visible','on'` ändern
- Datei: `Visualize_Terzband_Filter.m`
- Zeile: 74 - `'Visible','off'` auf `'Visible','on'` ändern

### Excel-Schreibfehler
Falls Excel-Export fehlschlägt:
- Prüfe Schreibrechte im `Excel/` Ordner
- Schließe geöffnete Excel-Dateien
- MATLAB benötigt Excel-Unterstützung (Windows) oder alternative Writer

### Fehlende Dateien
Das Skript überspringt fehlende Positionen automatisch und gibt Warnungen aus:
```
WARNUNG: Datei nicht gefunden für Variante_X, Position Y
```

---

## Literatur

- **IEC 61260:** Electroacoustics - Octave-band and fractional-octave-band filters
- **Lundeby et al. (1995):** "Uncertainties of Measurements in Room Acoustics"
- **Schroeder (1965):** "New Method of Measuring Reverberation Time"

---

## Autor

Erstellt für die Auswertung von Ultraschall-Raumimpulsantworten in Stuttgart.

**Letzte Aktualisierung:** 2026-01-06

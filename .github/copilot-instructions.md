# Copilot Instructions for PerformanceEngine

## Big Picture
- MATLAB workflows for combustion experiments: ingest Excel metadata, map to raw TDMS/TIF files, build a table object, then run pressure/photodiode/image processing and generate videos/statistics. Main pipeline scripts live at the repo root; CH4DDI/Engine_Code/New Code holds a newer engine-focused variant.
- Processing relies on a saved MATLAB object (`obj`) that wraps the data table and condition matrix. Configuration files define absolute data locations and must be updated per environment before running anything.

## Primary Pipeline (root)
- Update paths in [JournalFigures/Configuration.m](JournalFigures/Configuration.m) to point `DataSheetDirectory`, `DataTopfolder`, `DataObjSaveDir`, outputs, and `FieldNameExtended`. The active config at the bottom targets Methanol dual-fuel data and adds `Library_Matlab` + `tdms_package` to the path.
- Typical run order: (1) rename files if needed via [Step00_RenameFiles.m](Step00_RenameFiles.m); (2) build table and map raw files with [Step01_GenerateTableObject.m](Step01_GenerateTableObject.m); (3) process pressure/photodiode/HRR with [Step02_DataProcess.m](Step02_DataProcess.m); (4) compose condition videos via [Step03_MatrixVideo.m](Step03_MatrixVideo.m); (5) compute statistics with [Step04_StatisticalAnalysis.m](Step04_StatisticalAnalysis.m); (6) penetration boundary extraction in [Step05_Penetration.m](Step05_Penetration.m); (7) radial integration / flame metrics in [Step06_RadialIntegration.m](Step06_RadialIntegration.m).
- Scripts expect `Configuration()` on the MATLAB path; each script loads `Configs.DataObjReadDir` then writes back to `Configs.DataObjSaveDir`, so keep these in sync.

## Table and Condition Handling
- [Class_TableManagement.m](Class_TableManagement.m) reads Excel (Sheet1 data, Sheet2 conditions) with `ConditionMatchColumn` to align runs to conditions and builds run lists per condition. `Discard` column drives inclusion; `DataMatrix_runsNotMeetConditionAfterDiscard` holds mismatches.
- Use `ExtendTableWithDouble/ExtendTableWithChar` to add storage columns, and `LocateRawFileDirectory` to derive `File_Image`, `File_Pressure`, `File_Log`, `File_Photodiode` from date/run patterns under `Configs.DataTopfolder`.
- Condition tags for plotting/video come from `Class_TableManagement.ConditionTagGenerator`; boundary between table/condition matrices is critical when adding new derived metrics.

## Signal/Image Processing Patterns
- Pressure: [Library_Matlab/CLASS_PressureProcess.m](Library_Matlab/CLASS_PressureProcess.m) reads TDMS, applies hard-coded sampling rate, injection delay, voltage-to-bar scaling, computes corrected pressure, HRR, ignition delay, and optional plots. Adjust `SamplingRate`, `InjectionDelay`, or filtering routines there instead of per-script hacks.
- Photodiode: [Library_Matlab/CLASS_PhotodiodeProcess.m](Library_Matlab/CLASS_PhotodiodeProcess.m) expects `Gain_PD` and `InjectionDelay`; populated in [Step02_DataProcess.m](Step02_DataProcess.m) via `Table_PhotodiodeProcess`.
- Imaging: [Library_Matlab/CLASS_MultiPageTif.m](Library_Matlab/CLASS_MultiPageTif.m) loads stacks; [Step03_MatrixVideo.m](Step03_MatrixVideo.m) crops by `Nozzle_X`, resizes, overlays time/frame labels, and uses `CLASS_AxesHandleStore` for grids. Boundary/penetration extraction in [Step05_Penetration.m](Step05_Penetration.m) and radial integration in [Step06_RadialIntegration.m](Step06_RadialIntegration.m) assume cropped images and fixed `MaxLengthPixel`/`EdgeOffset`.
- HRR and ignition delay thresholds vary by fuel/type in `Table_HRRIgnitionDelay` inside [Step02_DataProcess.m](Step02_DataProcess.m); maintain these tables when adding new combustion modes.

## Outputs and Persistence
- `obj` is saved/loaded as a `.mat` file; each step mutates and re-saves it. Keep backward compatibility with existing `FieldNameExtended` fields to avoid load-time errors.
- Common figure/video/mat output directories are defined in `Configuration`; scripts rely on those absolute paths and will warn/fail if the folders are missing.

## CH4DDI Engine Variant (New Code)
- Configuration for engine TDMS batches is in [CH4DDI/Engine_Code/New Code/Configuration.m](CH4DDI/Engine_Code/New Code/Configuration.m); update Excel/data/output paths and ensure `Library_Matlab` is on the MATLAB path.
- [CH4DDI/Engine_Code/New Code/Library_Matlab/CLASS_TableManagement_Engine.m](CH4DDI/Engine_Code/New Code/Library_Matlab/CLASS_TableManagement_Engine.m) is a simplified table manager (no condition sheet) that locates multiple Takes per Set (`_Set##_*_TakeNN.tdms`) and stores them in `File_TDMS_AllTakes`.
- Engine workflow: run [Step01_GetTableObject.m](CH4DDI/Engine_Code/New Code/Step01_GetTableObject.m) to read Excel, extend fields, and locate TDMS files; run [Step02_PressureProcessEngine.m](CH4DDI/Engine_Code/New Code/Step02_PressureProcessEngine.m) to iterate sets and feed takes into `CLASS_PressureProcessEngine` (TODO: extend to derive corrected pressure/HRR). Diagnostics live in [Test_setup.m](CH4DDI/Engine_Code/New Code/Test_setup.m) and [Test_step2.m](CH4DDI/Engine_Code/New Code/Test_step2.m).

## Conventions and Gotchas
- Many scripts start with `clearvars; close all; clc;` and contain `keyboard` breakpoints—remove or continue execution manually when automating.
- Absolute Windows paths are hard-coded; adjust before committing or parameterize via `Configuration` to avoid environment-specific diffs.
- `tdms_package/Version_2p5_Final` is required for TDMS I/O; ensure it is on the path before calling any processing class.
- Adding new metrics: extend `FieldNameExtended`, update table extension in the relevant Step scripts, and persist via `obj.SetDataMatrix` to keep saves/loaders aligned.

## Engine Data Processing: TDMS Legacy Compatibility Rules (2026-01-05)
**Context:** CLASS_PressureProcessEngine must exactly match TDMS_Processing.m legacy behavior. See [VALIDATION_REPORT.md](CH4DDI/Engine_Code/New Code/VALIDATION_REPORT.md) for detailed analysis.

### Filtering (CRITICAL: Exact Match Required)
- **Pressure per-cycle:** 2nd-order Butterworth, `fc = 10000 Hz`, normalized as `fc_norm = fc/(7200*(RPM/60/2))`, applied with `filtfilt(b,a,P_raw)` (zero-phase).
  - Origin: TDMS_Processing.m lines 190-193
  - Implementation: CLASS_PressureProcessEngine.m, `PressureFilter` method
- **Pmean smoothing:** Savitzky-Golay (order=7, window=91) applied AFTER motoring correction, BEFORE HRR calculation.
  - Origin: TDMS_Processing.m line 326
  - **STATUS: MISSING in CLASS** → Add `P_corrected = sgolayfilt(P_corrected, 7, 91);`
- **aHRR smoothing:** Savitzky-Golay (order=7, window=91) applied AFTER forward-difference HRR calculation.
  - Origin: TDMS_Processing.m line 367
  - **STATUS: MISSING in CLASS** → Add `aHRR = sgolayfilt(aHRR, 7, 91);`
- **CA encoder smoothing:** Moving average, `b = ones(filter_width,1)/filter_width`, `filter_width = round((60/rpm/1800)/dt*5)`.
  - Origin: TDMS_Processing.m lines 168-171
  - Implementation: CLASS_PressureProcessEngine.m line 383

**Rule:** Do NOT modify filter parameters. Comment each with TDMS origin line number for traceability.

### Injection Signal Handling (CRITICAL: Di3 is Primary)
- **Channel priority:** Di3 (pilot diesel injection) is the earliest physical timing reference.
  - Primary: `Di3` (digital channel for pilot injection)
  - Fallback: `Untitled_3` (numerical alias for Di3)
  - Last resort: `Untitled_5` (main fuel injection, occurs later)
- **Current bug:** CLASS checks `Untitled_5` before `Di3` (wrong priority).
  - Fix: Reorder to check `Di3` → `Untitled_3` → `Untitled_5`
  - Log which channel is used for debugging

**Rule:** Always detect injection timing before AHRR calculation. Never allow AHRR features to precede injection in plots or metrics.

### Injection-Based AHRR Guard (NEW: Prevent Artifacts)
- **Physical principle:** Combustion cannot precede injection. Any pre-injection AHRR rise is a processing artifact (filtering edge effects, baseline drift).
- **Implementation:** For each cycle:
  1. Detect injection edge from Di3 (earliest rising edge in cycle)
  2. Map injection sample index to CA index using CA_Enc
  3. Zero all aHRR values before injection CA index
  4. Baseline-correct by subtracting pre-injection mean aHRR
- **Threshold for warnings:** Pre-injection aHRR > 10% of peak aHRR in that cycle
- **Origin:** Not in TDMS legacy (TDMS used injection-centered extraction which inherently prevented this)

**Rule:** Run injection guard BEFORE computing ensemble-averaged aHRR. Do not skip cycles—correct them via baseline subtraction.

### Per-Cycle Computation (CRITICAL: Compute Before Averaging)
- **Current bug:** CLASS computes CA10/CA50/CA90 from ensemble-averaged cumHRR (incorrect).
- **Correct method:**
  1. Loop over each cycle
  2. Compute aHRR, cumHRR, CA10, CA50, CA90 per cycle
  3. Apply SG filter to each cycle's aHRR
  4. Store arrays: `CA10_cycles(nCycles)`, `CA50_cycles(nCycles)`, `CA90_cycles(nCycles)`
  5. Compute mean and std: `CA10 = mean(CA10_cycles)`, `CA10_std = std(CA10_cycles)`
- **TDMS behavior:** Averaged pressure per-cycle before HRR computation (similar effect, proper variance handling).

**Rule:** Never compute metrics from already-averaged signals. Store per-cycle outputs in structured arrays. Perform averaging only at final reporting stage.

### CA Assignment Integrity (Encoder Quality Checks)
- **Encoder issues:** Noise can create duplicate CA samples or non-monotonic ordering.
- **Detection:** After extracting CA_cycle from encoder, run:
  ```matlab
  [CA_unique, uniq_idx] = unique(CA_cycle, 'stable');
  n_duplicates = numel(CA_cycle) - numel(CA_unique);
  if n_duplicates > 0
      warning('Cycle %d: %d CA duplicates (%.1f%%)', cc, n_duplicates, 100*n_dup/n_total);
  end
  ```
- **Correction:** Use `unique(..., 'stable')` to preserve temporal order. If >1% duplicates, warn (encoder quality issue).
- **Verification:** Assert `issorted(CA_unique)` and `all(diff(CA_unique) > 0)` (strictly monotonic).

**Rule:** Add CA uniqueness assertion per cycle. Never assume CA vectors are clean. Prefer interpolation onto reference CA grid if encoder data is poor.

### Z-Pulse and TDC Detection (Engine-Specific)
- **Z-pulse edge:** Use FALLING edge (1→0 transition) as TDC marker (confirmed empirically, +48 CAD offset applied).
- **TDC_shift parameter:** Adjusts cycle alignment in CA degrees (current: -68 CAD for engine tests).
- **Encoder channels:** Di7/Untitled_7 (Z-pulse), Di6/Untitled_6 (0.2° encoder pulses).

**Rule:** Z-pulse polarity is experiment-specific. Do not change without physical verification. The +48 CAD offset compensates for Z-pulse disk mechanical offset relative to thermodynamic TDC.

### Documentation Style
- **Imperative comments:** "Ensure...", "Do not...", "Apply...", "Zero..."
- **Physical justification:** Explain WHY (e.g., "Combustion cannot precede injection") not just WHAT.
- **TDMS traceability:** Reference legacy file and line number (e.g., `% TDMS_Processing.m:326`).
- **Warnings over silent fixes:** Log diagnostics when correction thresholds are exceeded.

**Example:**
```matlab
% === INJECTION-BASED AHRR GUARD (Validation Task 2) ===
% Physical principle: Combustion cannot precede injection
% TDMS legacy: Injection-centered extraction prevented this implicitly
% Implementation: Zero aHRR before Di3 injection edge per cycle
```

## Quick Start for Agents
- Edit `Configuration.m` for your dataset paths and add required folders to the MATLAB path.
- Run the Step scripts sequentially, reloading `obj` each time; verify intermediate tables/outputs before progressing.
- When debugging data issues, inspect the table object (`obj.DataMatrix`, `obj.ConditionMatrix`) to confirm file paths, discard flags, and run-number alignment.
- **For engine data:** See [VALIDATION_REPORT.md](CH4DDI/Engine_Code/New Code/VALIDATION_REPORT.md) for TDMS compatibility rules and validation tests.

Feedback welcome—call out unclear sections so we can refine these notes.

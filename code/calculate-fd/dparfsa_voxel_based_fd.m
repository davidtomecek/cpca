% =========================================================================
% DPARSFA HEADLESS PIPELINE: VOXEL-WISE & REGION-BASED MOTION ANALYSIS
% =========================================================================
% Methodology: Yan et al. (2013) NeuroImage.
% Implemented via programmatic DPABI/DPARSF API.

clear; clc;

% -------------------------------------------------------------------------
% 1. PATHS & DIRECTORY SETUP
% -------------------------------------------------------------------------
% MODIFY THESE PATHS to match your local installation or server layout
SPM_Path   = '/hydra/hydra_io/vypocty/tomecek/toolbox/spm12';
DPABI_Path = '/hydra/hydra_io/vypocty/tomecek/toolbox/DPABI_V9.0_250415';

% Add frameworks dynamically to the top of your MATLAB search path
if ~exist('spm', 'file'), addpath(genpath(SPM_Path)); end
if ~exist('DPARSFA_run', 'file'), addpath(genpath(DPABI_Path)); end

% Initialize core configuration structure
Cfg = [];

% Set the base directory containing your 'FunImg' and 'T1Img' folders
Cfg.DataProcessDir = '/hydra/hydra_io/vypocty/tomecek/cpca/data/TEMP'; 

% -------------------------------------------------------------------------
% 2. AUTOMATED SUBJECT EXTRACTION (CONVENTION-OVER-CONFIGURATION)
% -------------------------------------------------------------------------
Cfg.StartingDirName   = 'FunImg';      % Raw 4D functional directory
Cfg.StructuralDirName = 'T1Img';       % Raw high-res structural directory
Cfg.SubjectIDPattern  = 'ESO_*';       % String wildcards matching folders

% Automatically populate subject arrays based on directory layout
SubjectDirs = dir(fullfile(Cfg.DataProcessDir, Cfg.StartingDirName, Cfg.SubjectIDPattern));
if isempty(SubjectDirs)
    error('Pipeline Aborted: No subject folders matching pattern "%s" found in %s/FunImg/', ...
        Cfg.SubjectIDPattern, Cfg.DataProcessDir);
end
Cfg.SubjectID = {SubjectDirs.name}';

% -------------------------------------------------------------------------
% 3. FUNCTIONAL PARAMETERS & ABSOLUTE SLICE TIMING (BIDS-STYLE)
% -------------------------------------------------------------------------
Cfg.TR = 2.0;                          % Repetition Time (in seconds)
Cfg.TimePoints = 0;                    % 0 forces auto-detection of 4D volume length

Cfg.IsSliceTiming = 1;                 % Activate slice-timing correction
Cfg.SliceTiming.SliceNumber = 37;      % Total physical slice count

% PASTE YOUR REAL SLICE ACQUISITION TIMES BELOW (Units MUST match your TR scale)
% Decimals force the underlying SPM engine into absolute temporal interpolation.
% Length of array MUST perfectly match Cfg.SliceTiming.SliceNumber.
Cfg.SliceTiming.SliceOrder = [
                0,
                1.015,
                0.0525,
                1.07,
                0.1075,
                1.1225,
                0.16,
                1.1775,
                0.2125,
                1.23,
                0.2675,
                1.2825,
                0.32,
                1.3375,
                0.375,
                1.39,
                0.4275,
                1.445,
                0.48,
                1.4975,
                0.535,
                1.5525,
                0.5875,
                1.605,
                0.6425,
                1.6575,
                0.695,
                1.7125,
                0.7475,
                1.765,
                0.8025,
                1.82,
                0.855,
                1.8725,
                0.91,
                1.925,
                0.9625  ];

% CRITICAL: When SliceOrder is time, ReferenceSlice MUST be an absolute time point.
% 0.0 anchors temporal alignment to the exact start of the volume capture window.
Cfg.SliceTiming.ReferenceSlice = 0.0; 

% -------------------------------------------------------------------------
% 4. PREPROCESSING BLOCKS (RIGID METRIC SETUP)
% -------------------------------------------------------------------------
Cfg.IsRealign = 1;                     % Calculate rigid-body translations/rotations

% Skip coregistration/segmentation for motion-only extraction (non-interactive)
Cfg.IsCoregister = 0;                  % Skip co-registration (avoids GUI dialog)
Cfg.IsSegment    = 0;                  % Skip segmentation (not needed for FD)

% -------------------------------------------------------------------------
% 5. ACTIVATE VOXEL-SPECIFIC HEAD MOTION (YAN ET AL. 2013 METHODOLOGY)
% -------------------------------------------------------------------------
% This creates spatial 4D voxel displacement arrays instead of abstract global estimates.
Cfg.IsCalVoxelSpecificHeadMotion = 1;  

% -------------------------------------------------------------------------
% 6. SAFEGUARDS: CONTROL OVERLAPPING STEP INTERACTION
% -------------------------------------------------------------------------
% CRITICAL: Spatial smoothing destroys physical voxel boundary coordinates.
% Keep at 0 to guarantee strict math execution of structural voxel arrays.
Cfg.IsSmooth = 0;                      
Cfg.IsFilter = 0;                      % High/Low pass filtering skipped for raw extraction

% -------------------------------------------------------------------------
% 7. COMPUTATIONAL RESOURCES & RUNTIME EXECUTION
% -------------------------------------------------------------------------
Cfg.IsParallel = 1;                    % Enable MATLAB Parallel Computing Toolbox
Cfg.ParallelWorkersNumber = 4;         % Number of allocated CPU cores/workers

% Announce configuration details to standard terminal output
fprintf('=======================================================\n');
fprintf('  LAUNCHING HEADLESS DPARSFA PIPELINE\n');
fprintf('  Processing Workspace: %s\n', Cfg.DataProcessDir);
fprintf('  Detected Subjects:    %d\n', length(Cfg.SubjectID));
fprintf('=======================================================\n');

% Send configuration properties directly to the processing loop execution engine
DPARSFA_run(Cfg);

fprintf('\nPipeline processing loop terminated successfully!\n');
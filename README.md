# EEG-Preprocessing-with-ADJUST

The main script is 'preprocess_signals.m'
Step 1: Ensure EEGLAB folder is available in Matlab paths or place the EEGLAB folder in the same repository. ADJUST plugin should be installed as well.
Step 2: Create a directory 'signals\raw\' and place EEG dataset in EEGLAB format (*.set) or change the path within the script
Step 3: Change preprocessing configurations (Optional)

interface_ADJ_mod.m is a modified version of interface_ADJ.m available in the ADJUST plugin. This version disable the interactive UI responsible of manual rejection of contaminated ICs.

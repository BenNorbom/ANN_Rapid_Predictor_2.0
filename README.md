# ANN Rapid Predictor 2.0

Computational models are powerful tools that can enable the optimization of deep brain stimulation (DBS). To enhance the clinical practicality of these models, their computational expense and required technical expertise must be minimized. An important aspect of DBS models is the prediction of neural activation in response to electrical stimulation. Existing rapid predictors of activation simplify implementation and reduce prediction runtime, but at the expense of accuracy. We sought to address this issue by leveraging the speed and generalization abilities of artificial neural networks (ANNs) to create a novel predictor of neural fiber activation in response to DBS.

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/BenNorbom/ANN_Rapid_Predictor_2.0.git
   cd ANN_Rapid_Predictor_2.0
   ```

2. Pull example tracks:
   To run the quick start examples, you'll need the tract files. You can pull just these files (~246 KB) without downloading the entire large file history:
   ```bash
   git lfs pull -I "example_tracks/**"
   ```

3. **(Important)** Pull large model and electrode files:
   This repository uses Git LFS for large files. Run the following command to download all remaining necessary files (approx. 669 MB):
   ```bash
   git lfs pull
   ```

   **Selective Download:**
   If you only need specific electrode files (to save bandwidth/disk space), you can pull them individually using the `-I` (include) flag. For example:
   ```bash
   git lfs pull -I "electrodes/directed/monopolar/bsc_directional_anisotropic_monopolar_01(a,b,c)2(-a,b,c)3.txt"
   ```

4. Create a Python 3.10 virtual environment and install dependencies:
   ```bash
   python -m venv .venv
   .venv\Scripts\activate      # Windows
   # source .venv/bin/activate  # Linux/Mac
   pip install -r requirements.txt
   ```

> **Note:** Python 3.10 is required for TensorFlow 2.15 compatibility.

## Quick Start

### 1. Run the ANN prediction

```bash
python run/dti_ann_LUT.py \
    electrodes/medtronic_3387/monopolar/3387_anisotropic_monopolar_01-23.txt \
    example_tracks/whole_brain.txt \
    models/ann_19_reg \
    run/results.json \
    ssd dti 38 50 30 reg
```

See [run/README.md](run/README.md) for detailed argument descriptions.

### 2. Visualize the results

```bash
python graphing/plot_tracts_fast.py \
    --tract example_tracks/whole_brain.txt \
    --results run/results.json \
    --output output_viz/ \
    --activation_threshold 3.0 \
    --electrode_center 38 50 30 \
    --electrode_config 01-23 \
    --downsample 3
```

See [graphing/README.md](graphing/README.md) for all visualization options.

### 3. Bundle-aware visualization (optional)

If you have specific tract bundles (e.g. DRTT, ML, PRT) alongside a whole-brain
file, you can merge them and visualize each bundle in a distinct colour:

```bash
# Merge whole-brain tracts with named bundles
python move_whole_brain_tracts.py \
    --source example_tracks/whole_brain.txt \
    --center 38 50 30 \
    --bundles L_DRTT_voxel.txt L_ML_voxel.txt L_PRT_voxel.txt \
    --output merged_tracts.txt

# Run prediction on the merged file
python run/dti_ann_LUT.py \
    electrodes/medtronic_3387/monopolar/3387_anisotropic_monopolar_01-23.txt \
    merged_tracts.txt models/ann_19_reg run/results.json \
    ssd dti 38 50 30 reg

# Visualize with bundle colouring
python graphing/plot_tracts_bundles.py \
    --tract merged_tracts.txt \
    --results run/results.json \
    --manifest merged_tracts_manifest.json \
    --output output_viz_whole_brain/ \
    --electrode_center 38 50 30 \
    --electrode_config 01-23
```

## Project Structure

```
ANN_Rapid_Predictor_2.0/
├── electrodes/          # Electrode voltage field data (FEM exports)
├── example_tracks/      # Example tractography files
├── graphing/            # Visualization scripts
│   ├── plot_tracts_fast.py      # Fast 3D visualization (100k+ fibers)
│   ├── plot_tracts.py           # Original visualization with per-fiber rendering
│   └── plot_tracts_bundles.py   # Bundle-aware colour-coded visualization
├── models/              # Pre-trained ANN models
├── run/                 # Prediction scripts
│   ├── dti_ann_LUT.py           # Main ANN prediction script
│   └── process_DTI.py           # Tract processing & spline interpolation
├── requirements.txt     # Python dependencies
└── README.md
```

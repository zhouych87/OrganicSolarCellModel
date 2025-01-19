markdown复制
# Python Scripts for Molecular Dynamics Analysis

This repository contains a collection of Python scripts designed to analyze molecular dynamics simulation data, specifically for the eC9 system. The scripts cover various stages of data processing, from feature extraction to clustering and visualization. Below is a detailed description of each script and its functionality.

## 1. ExtractStackingFeature.py

### Purpose
Extracts various features from molecular dynamics trajectories, including molecular pair stacking indices, segment average distances, maximum and minimum distances, stacking angles, segment and molecular centroid distances, and many internal molecular distance features.

### Input
- `topology_file`: "eC9.tpr" - A dynamics-generated .tpr file.
- `trajectory_file`: "eC9.trr" - A dynamics-generated .trr file.

### Output
- `eC9_AA.xlsx`, `eC9_AD.xlsx`, `eC9_DD.xlsx`, `satisfying_pairsAA`, `satisfying_pairsAD`, `satisfying_pairsDD` - Excel files containing the extracted features.

### Remarks
- Uses GROMACS commands to generate and convert files:
  - `gmx grompp -f opt.mdp -c confout.gro -p 3PM6eC9.top -o md_xtc.tpr` - Regenerate a suitable .tpr file.
  - `gmx trjconv -f traj.trr -s md_xtc.tpr -o opt.xtc -pbc mol` - Convert to .xtc format.

## 2. ClusterKmeans.py

### Purpose
Performs K-means clustering on selected features from the extracted data.

### Input
- `file_path`: 'eC9AD.xlsx' - An Excel file containing selected features (3 columns) from eC9_AD.xlsx.

### Output
- `eC9AD_Feature_Scaled.xlsx` - Scaled features.
- `eC9AD_Cluster.xlsx` - Clustering results without normalization.
- `eC9AD_Cluster_Scaled.xlsx` - Clustering results with normalization.
- `eC9AD_Cluster_fit.xlsx` - Data set close to the top-left category, with original unscaled features.
- `eC9AD_Cluster_Center.xlsx` - Information on the 5 closest and farthest points from each cluster center, with original unscaled features.
- `ClusterCenterInformation.xlsx` - Coordinates of each cluster center for graphical display.

### Remarks
- Uses GROMACS command to convert trajectory files:
  - `gmx trjconv -f traj.trr -s md_xtc.tpr -o opt.xtc -pbc mol`

## 3. ClusterKmeans3dplot.py

### Purpose
Generates 3D plots of the clustering results.

### Input
- `Cluster_Scaled`: 'eC9AD_Cluster_Scaled.xlsx' - Scaled clustering data.
- `ClusterCenterInformation`: 'ClusterCenterInformation.xlsx' - Cluster center coordinates.
- `Molecule`: 'eC9' - Molecule identifier.

### Output
- A 3D plot image without a legend.

### Remarks
- Automatically colors the plot from deep purple to yellow.
- Uses GROMACS command to convert trajectory files:
  - `gmx trjconv -f traj.trr -s md_xtc.tpr -o opt.xtc -pbc mol`

## 4. ClusterKmeans3dplot_legend_only.py

### Purpose
Generates a legend-only image for the 3D plots.

### Input
- `Cluster_Scaled`: 'eC9AD_Cluster_Scaled.xlsx' - Scaled clustering data.
- `ClusterCenterInformation`: 'ClusterCenterInformation.xlsx' - Cluster center coordinates.
- `Molecule`: 'eC9' - Molecule identifier.

### Output
- A legend image.

### Remarks
- Automatically colors the legend from deep purple to yellow.
- Uses GROMACS command to convert trajectory files:
  - `gmx trjconv -f traj.trr -s md_xtc.tpr -o opt.xtc -pbc mol`

## 5. Random.py

### Purpose
Generates random residue pair indices for further analysis.

### Input
- `file_path`: 'Random.xlsx' - A database of specific class residue pair indices.
- `num_samples`: 9 - Number of random indices to extract.
- `random_seed`: 8 - Random seed definition.

### Output
- `Grep.xlsx` - 0 residue pair indices (9 rows, 2 columns) for structure quantification and job submission file generation.
- `Random_Randomout.xlsx` - 9 random residue pair indices (9 rows, 4 columns) for comparison and effective distance calculation.

### Remarks
- Uses GROMACS command to convert trajectory files:
  - `gmx trjconv -f traj.trr -s md_xtc.tpr -o opt.xtc -pbc mol`

## 6. Grepgjf_batch.py

### Purpose
Batch generates .gjf files for quantum chemistry calculations.

### Input
- `tpr_file_path`: 'eC9.tpr' - Trajectory file.
- `trr_file_path`: 'eC9.trr' - Trajectory file.
- `df`: pd.read_excel('Grep.xlsx', header=None) - Residue pair information.

### Output
- .gjf files for all specified residue pairs in Grep.xlsx.

### Remarks
- Uses GROMACS command to convert trajectory files:
  - `gmx trjconv -f traj.trr -s md_xtc.tpr -o opt.xtc -pbc mol`

## 7. Calculate_Distance_Couple.py

### Purpose
Calculates effective distances for specified residue pairs.

### Input
- `random_out`: pd.read_excel('Randomout.xlsx', header=None) - Residue pair indices from Random_Randomout.xlsx.
- `aa_feature_scaled`: pd.read_excel('eC9AD_Feature_Scaled.xlsx') - Scaled features.

### Output
- `Scaled_Feature_Select.xlsx` - Normalized 3D data for all specified residue pairs.
- `Scaled_Feature_Select_Distance.xlsx` - Calculated effective distance data.

### Remarks
- Uses GROMACS command to convert trajectory files:
  - `gmx trjconv -f traj.trr -s md_xtc.tpr -o opt.xtc -pbc mol`

## Usage

1. **Feature Extraction**
   - Run `ExtractStackingFeature.py` to extract features from your molecular dynamics data.
   - Example: `python ExtractStackingFeature.py`

2. **Clustering**
   - Run `ClusterKmeans.py` to perform K-means clustering on the extracted features.
   - Example: `python ClusterKmeans.py`

3. **Visualization**
   - Run `ClusterKmeans3dplot.py` to generate 3D plots of the clustering results.
   - Run `ClusterKmeans3dplot_legend_only.py` to generate a legend for the plots.
   - Example: `python ClusterKmeans3dplot.py`

4. **Random Index Generation**
   - Run `Random.py` to generate random residue pair indices.
   - Example: `python Random.py`

5. **Batch .gjf Generation**
   - Run `Grepgjf_batch.py` to generate .gjf files for quantum chemistry calculations.
   - Example: `python Grepgjf_batch.py`

6. **Distance Calculation**
   - Run `Calculate_Distance_Couple.py` to calculate effective distances for specified residue pairs.
   - Example: `python Calculate_Distance_Couple.py`

## Dependencies

- Python 3.x
- GROMACS
- Pandas
- NumPy
- Matplotlib (for plotting)

## Contributing

Contributions are welcome! Please open an issue or submit a pull request if you have any suggestions or improvements.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

For any questions or further information, please contact [zlzhang@gzu.edu.cn].

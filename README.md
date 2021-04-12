# BM2-MD
BM2 collaboration - Simulation analysis scripts

# Overview
All scripts take two inputs: a gromacs simulation input file [TPR] and a gromacs trajectory [XTC/TRR]. Additional details on what each script calcualtes are provided below. For help, execute a script with the `-h` flag.
When projections are used, they are based on the Ca coordinates of the relevant residues.

# Citations
MDAnalysis is used in these scripts. References for MDAnalysis are:

Gowers et al. S. Benthall and S. Rostrup ed. Proceedings of the 15th Python in Science Conference, Austin, TX, 2016: 98.

Michaud-Agrawal et al. J. Comput. Chem. 2011(32): 2319.

MDAnalysis.analysis.hbonds.HydrogenBondAnalysis is used in water_hydrogen_bond_network.py. The reference (for a later version) is:
Smith et al. Phys. Chem. Chem. Phys. 2019(21): 9845.

MDAnalysis.analysis.rms.RMSD is used in protein_RMSD.py. The references are:
Theobald Acta Crystallographica A 2005(61): 478.
Liu et al. J. Comput. Chem. 2010(31): 1561.

# water_structure.py (Fig. 2d-e, Fig S3a, Fig S7, Fig S8)
Example: ./water_structure.py md.tpr md.xtc
Calculates the water density (Fig. 2d-e), OH order parameter (Fig 5e), HH order parameter, and cosine(theta_{OH}) distribution (Fig 5d, S7, S8). Can be modified to calculate these values as a function of time (by binning in time), or to calculate the properties over the entire channel. By default, the output has each property as a function of the channel-axis coordinate, Z, with bin-centers at -32:4:32.

# water_transport.py (Fig 3g)
Example: ./water_transport.py md.tpr md.xtc
Returns the number of waters that are transported through the channel as a function of time (a water is recorded at the time where it exits the channel). The transported.txt output contains a column for transport from the C to N terminal end, and a column for transport from the N to C terminal end (Fig 3g). The script also returns the passage time for each water molecule, which is the time a water molecule spent in the channel during transport (mentioned in the manuscript text).

# water_MSD.py (Fig S3b-c)
Example: ./water_MSD.py md.tpr md.xtc
Returns the mean-squared displacement of water molecules in the channel, as a function of channel-axis coordinate (-32:4:32 range). By default, MSD is calculated for time 0:10ps:200ps. This can be adjusted as needed. Each water molecules' contrubition is added based on the bin where the water molecule residues at \tau=0 (i.e., water molecules can move between bins and will still be included in the analysis).

# water_Crot.py (Fig 3e-f, S3a)
Example: ./water_Crot.py md.tpr md.xtc
Returns 0.5*<3*\mu(t)\cdot\mu(t+\tau) - 1>_t and <\mu(t)\cdot\mu(t+\tau)>_t for the dipole vector (\mu) of water molecules in the channel, as a function of channel-axis coordinate (-32:4:32 range). By default, calculated for \tau=0:10ps:200ps. This can be adjusted as needed. Each water molecules' contribution is added based on the bin where the water molecule resides at \tau=0 (i.e., water molecules can move between bins and will still be included in the analysis).

# water_hydrogen_bond_network.py (Fig. 7b-d)
Example:
./water_hydrogen_bond_network.py \\  
--directory_for_output /network_analysis_output_directory \\  
--title_of_job network_analysis_job \\  
--dictionary_of_simulations simulations_dictionary.json \\  
--selected_simulations P01 P4 \\  
--dictionary_of_colors colors_dictionary.json \\  
--minimum_axis_value -20 \\  
--maximum_axis_value 22 \\  
--axis_increments 2 \\  
--frame_stride 10 \\  
--minimum_y_for_plot -21 \\  
--maximum_y_for_plot 21 \\  
--residue_numbers_for_label 8 12 16 19 23 27 &  
Returns statistics on the hydrogen bonds, hydrogen bond directionality, and bottlenecks between adjacent channel axis slices. For this script, the slices are centered at -20, -18, ...22, with each bin having width 2. This is updated from the script used for the Communications Biology figure, where the first and last slice had length 1, with all other slices having length 2 (so this only affects results at the edges). A printout and plots are generated, as well as a .npy data file. The plots were then processed in Adobe Illustrator for clarity and standardization (e.g., setting projections to be standardized on all plots, to a higher-strided projection value used earlier). Note all frames analyzed are included in all analyses (i.e., frames where there are no waters in a slice are not excluded- preliminary checks indicate such exclusion does not drastically alter results). There is a more recent version of the hydrogen bond module (MDAnalysis.analysis.hydrogenbonds.hbond_analysis) which, to the best of our knowledge, did not exist when this analysis was developed. Updating to the newer module may be of interest in the future. One observation that was made after this workflow was developed is that the standard settings for this module appear to be overly generous, giving more hydrogen bonds than expected for protein atoms when analyzing protein-water hydrogen bonds (https://docs.mdanalysis.org/stable/documentation_pages/analysis/hbond_analysis.html explains the angle cutoff may be why). Running for a different angle cutoff indicated qualitative trends did not change drastically. While revising the cutoffs may be of interest in the future, the default settings' being generous here may be helpful, because with striding, this helps ensure hydrogen bonds are not missed. The argparse section (obtained by running "./water_hydrogen_bond_network.py -h") contains more information on the different arguments.

# protein_minimum_distance.py (Fig.S1b)
Example:
./protein_minimum_distance.py \\  
--selected_residue_numbers 8 12 16 19 23 27 --residue_numbers_for_label 8 12 16 19 23 27 \\  
--directory_for_output /minimum_distance_output_directory \\  
--title_of_job minimum_distance_job \\  
--dictionary_of_simulations simulations_dictionary.json \\  
--selected_simulations P01 P4 \\  
--dictionary_of_colors colors_dictionary.json \\  
--frame_stride 1 \\  
--minimum_x_for_plot 0 \\  
--maximum_x_for_plot 18 \\  
--minimum_y_for_plot -18 \\  
--maximum_y_for_plot 18 \\  
--dictionary_of_structures structures_dictionary.json &  
Returns minimum distances between heavy atoms of the same residue on diagonally placed side chains. A printout and plot are generated, as well as a .json data file. The plot was then processed in Adobe Illustrator for clarity and standardization (e.g., setting projections to be standardized on all plots, to a higher-strided projection value used earlier). The argparse section (obtained by running "./protein_minimum_distance.py -h") contains more information on the different arguments.

# protein_RMSD.py (Fig. S1c-d)
Example:
./protein_RMSD.py \\  
--directory_for_output /RMSD_output_directory \\  
--title_of_job RMSD_job \\  
--dictionary_of_simulations simulations_dictionary.json \\  
--dictionary_of_colors colors_dictionary.json \\  
--dictionary_of_selections selections_dictionary.json \\  
--frame_stride 10 \\  
--minimum_y_for_plot 0.0 \\  
--maximum_y_for_plot 5.0 \\  
--minimum_x_for_plot 0.0 \\  
--maximum_x_for_plot 135.0 &  
Returns RMSDs of selected portions of the protein over the trajectory, each relative to a reference. A printout and plots are generated, as well as a .json data file. The plots were then processed in Adobe Illustrator for clarity and standardization (e.g., changing colors). The argparse section (obtained by running "./protein_RMSD.py -h") contains more information on the different arguments.



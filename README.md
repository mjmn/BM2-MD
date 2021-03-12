# BM2-MD
BM2 collaboration - Simulation analysis scripts

# Overview
All scripts take two inputs: a gromacs simulation input file [TPR] and a gromacs trajectory [XTC/TRR]. Additional details on what each script calcualtes are provided below. For help, execute a script with the `-h` flag.

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
./water_hydrogen_bond_network.py \
--directory_for_output /output_directory \
--title_of_job network_analysis_job \
--dictionary_of_simulations simulations_dictionary.json \
--selected_simulations P01 P4 \
--dictionary_of_colors colors_dictionary.json \
--minimum_axis_value -20 \
--maximum_axis_value 22 \
--axis_increments 2 \
--frame_stride 10 \
--minimum_y_for_plot -21 \
--maximum_y_for_plot 21 \
--residue_numbers_for_label 8 12 16 19 23 27 &
\nReturns statistics on the hydrogen bonds, hydrogen bond directionality, and bottleneecks between adjacent channel axis slices (-20:2:22 range, with first and last bins having length of 1, which may be useful for axis scaling so bulk data does not cause scaling issues and may be problematic in other cases). A printout and plots are generated. The plots were processsed in Adobe Illustrator for clarity. Note all frames are included in all analyses, and also that there is a more recent version of the hydrogen bond module (MDAnalysis.analysis.hydrogenbonds.hbond_analysis) which, to the best of our knowledge, did not exist when this analysis was developed


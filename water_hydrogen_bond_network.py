import argparse
import glob
from itertools import groupby
import json
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis
import MDAnalysis.analysis.hbonds
import numpy as np
import numpy.linalg as npla
import os
import scipy
from scipy import stats
import warnings

#############
#For each replicate of the simulation trajectories specified by simulation_selection_name, calculate the H bonds between channel axis slices, the H bond bottlenecks at each channel axis slice, and the directionality metric at each channel axis slice, as well as bottlenecks across the channel
#############
def network_analysis(simulation_selection_name, information_dictionary, args, ch_ax_list, analysis_output):

    #############
    #Obtain the individual replicate directories for each directory (https://stackoverflow.com/questions/14506491/php-glob-pattern-match-for-arbitray-number-of-digits) and document in output file created
    #############
    print(f"On {simulation_selection_name}") #Use f-strings to format (https://realpython.com/python-f-strings/)
    analysis_output.write(f"=================\n{simulation_selection_name}\n")
    replicate_directories_1_digit = glob.glob(f"{information_dictionary[simulation_selection_name]['parent_directory']}[0-9]")
    replicate_directories_2_digit = glob.glob(f"{information_dictionary[simulation_selection_name]['parent_directory']}[0-9][0-9]") #This is written for a directory system where there is a numerical suffix for each replicate, and there are up to 100 replicates
    replicate_directories = sorted(replicate_directories_1_digit) + sorted(replicate_directories_2_digit)
    analysis_output.write(f"*************\nAnalysis for {simulation_selection_name}\nParent directory {information_dictionary[simulation_selection_name]['parent_directory']}\n")
    analysis_output.write("Replicate directories are\n")
    for rep_d in replicate_directories:
        analysis_output.write(f"{rep_d}\n")

    #############
    #Set up data structures to hold simulation data analyzed for each replicate- H bond metrics, frame count, and projections
    #############
    #np array where each entry is bottleneck count for one frame analyzed
    overall_bottleneck_counts = np.empty(0)

    #np array where each entry is averaged bottleneck count for one replicate's frames analyzed
    av_overall_bottleneck_counts = np.empty(0)

    #Dictionary where keys are channel axis coordinates and values are the average count for each replicate of hydrogen bonds between water molecules in this channel axis coordinate bin and those in the adjacent channel axis coordinate bin closer to the C-terminus
    #For this and the two dictionaries below, set up the values for each channel axis coordinate key to be an empty np array. Each np array will be populated with replicate-level averages from the np arrays of data that replicate_network_analysis will generate.
    av_overall_hb_dictionary = dict.fromkeys(ch_ax_list[:-1], np.empty(0)) #https://www.programiz.com/python-programming/methods/dictionary/fromkeys

    #Dictionary where keys are channel axis coordinates and values are the average for each replicate of the hydrogen bond directionality metric between water molecules in this channel axis coordinate bin and those in the adjacent channel axis coordinate bin closer to the C-terminus
    av_overall_hb_diff_dictionary = dict.fromkeys(ch_ax_list[:-1], np.empty(0))

    #Dictionary where keys are channel axis coordinates and values are the percentage for each replicate of frames in which there is a bottleneck (no hydrogen bonds present) between water molecules in this channel axis coordinate bin and those in the adjacent channel axis coordinate bin closer to the C-terminus
    av_overall_bottleneck_dictionary = dict.fromkeys(ch_ax_list[:-1], np.empty(0))

    #If the user wants to label residues on the plot, set up a dictionary where keys are these residues' numbers and values are np arrays of channel axis projections
    if (args.residue_numbers_for_label != None):
        ch_ax_proj_dictionary = dict.fromkeys(args.residue_numbers_for_label, np.empty(0))

    #Set up frame counting
    total_num_fr = 0 #Total frame count, to return and compare to frame count for other simulation conditions analyzed
    frame_count_array = np.empty(0) #np array that will hold lengths of each replicate's data array, so a warning can be printed if they differ

    #############
    #For each replicate, call replicate_network_analysis to carry out analysis. Record the results in the data structures collecting data for each replicate
    #############
    for replicate_directory in replicate_directories:
        print(f"On {replicate_directory}")
        if (os.path.isfile(f"{replicate_directory}/md_1.xtc")):

            #############
            #Call replicate_network_analysis to run analysis for this replicate
            #############
            bottleneck_counts, hb_dictionary, hb_diff_dictionary, bottleneck_dictionary, num_fr, res_proj_dictionary, num_to_name_dictionary = replicate_network_analysis(args, ch_ax_list, replicate_directory, analysis_output)
        
            #############
            #Update the overall data structures
            #############
            #Add np array of all bottleneck counts to growing overall_bottleneck_counts np array
            overall_bottleneck_counts = np.append(overall_bottleneck_counts, bottleneck_counts)

            #Add average bottleneck count for this replicate to the growing av_overall_bottleneck_counts np array
            av_overall_bottleneck_counts = np.append(av_overall_bottleneck_counts, np.average(bottleneck_counts))

            #For each channel axis slice, add the average hydrogen bond count, average directionality metric, and bottleneck percentage for this replicate to the np array for this channel axis slice in the dictionaries av_overall_hb_dictionary, av_overall_hb_diff_dictionary, and av_overall_bottleneck_dictionary, respectively.
            for k_ch in av_overall_hb_dictionary.keys():
                av_overall_hb_dictionary[k_ch] = np.append(av_overall_hb_dictionary[k_ch], np.average(hb_dictionary[k_ch]))
                av_overall_hb_diff_dictionary[k_ch] = np.append(av_overall_hb_diff_dictionary[k_ch], np.average(hb_diff_dictionary[k_ch]))
                av_overall_bottleneck_dictionary[k_ch] = np.append(av_overall_bottleneck_dictionary[k_ch], np.average(bottleneck_dictionary[k_ch])) #Averaging the 1s/0s gives the percentage of frames with a bottleneck

            #Record frame count for this replicate, adding to total_num_fr and entering the value into frame_count_array
            total_num_fr = total_num_fr + num_fr
            frame_count_array = np.append(frame_count_array, num_fr)

            #If the user wants to label residues on the plot, add average channel axis coordinate projection for each residue of interest for this replicate (averaging the data stored in res_proj_dictionary) to the np array in ch_ax_proj_dictionary corresponding to this residue that holds the average projection for each replicate
            if (args.residue_numbers_for_label != None):
                for residue_to_label_in_plot in args.residue_numbers_for_label:
                    ch_ax_proj_dictionary[residue_to_label_in_plot] = np.append(ch_ax_proj_dictionary[residue_to_label_in_plot], np.average(res_proj_dictionary[residue_to_label_in_plot]))

    #############
    #Compute and report summary statistics across replicates
    #############
    #After all replicates have been analyzed, check that all of them have the same count of observations. If they do not, print out a warning that the average across frames will differ from the average across replicates.
    number_of_distinct_obs_counts = len(np.unique(frame_count_array)) #Obtain count of unique values (https://kite.com/python/answers/how-to-count-frequency-of-unique-values-in-a-numpy-array-in-python)
    if (number_of_distinct_obs_counts > 1):
         warnings.warn("***Please note observation counts are not the same for all replicates. Be aware that the cross-replicate average, which is reported, will be different from the average obtained by weighting each frame equally.")
         analysis_output.write("***Please note observation counts are not the same for all replicates. Be aware that the cross-replicate average, which is reported, will be different from the average obtained by weighting each frame equally.\n")
         print(frame_count_array)

    #Calculate the average count of bottlenecks across frames- average np array of average bottleneck counts, this is an average of the replicate-level averages
    average_bott = np.average(av_overall_bottleneck_counts)
    analysis_output.write(f"AVERAGE {simulation_selection_name} ACROSS REP BOTTLENECK PER FR COUNT {average_bott:2f}\n")

    #Calculate the cross-replicate average H bond count, directionality metric, and bottleneck percentage for each channel axis slice
    for ch_a_k in ch_ax_list[:-1]:
        analysis_output.write(f"Axis {ch_a_k}: average hb count {np.average(av_overall_hb_dictionary[ch_a_k]):3f}, average directionality metric {np.average(av_overall_hb_diff_dictionary[ch_a_k]):3f}, average bottleneck percent {np.average(av_overall_bottleneck_dictionary[ch_a_k]):3f}\n\n")
    analysis_output.write("\n\n")

    #If projections are wanted, calculate the average projection of residues onto the channel axis
    if (args.residue_numbers_for_label != None):
        for r in ch_ax_proj_dictionary.keys():
            analysis_output.write(f"Residue {num_to_name_dictionary[r]}: average projection {np.average(ch_ax_proj_dictionary[r]):2f}, sem {scipy.stats.sem(ch_ax_proj_dictionary[r]):2f}\n")

    #If projections are not wanted, set ch_ax_proj_dictionary to be None
    if (args.residue_numbers_for_label == None):
        ch_ax_proj_dictionary = None

    analysis_output.write("\n\n")

    #Return bottleneck counts np array, frame count value, residue number-to-label mapping dictionary, channel axis projection dictionary, and dictionaries with each channel axis slice's H bond count, H bond directionality metric, and bottleneck percentage
    return overall_bottleneck_counts, total_num_fr, num_to_name_dictionary, ch_ax_proj_dictionary, av_overall_hb_dictionary, av_overall_hb_diff_dictionary, av_overall_bottleneck_dictionary

#############
#Given a list of residue numbers of water molecules to include (wat_res_num_list), returns a MDAnalysis selection
#byres is used to obtain entire water molecules (https://www.mdanalysis.org/docs/documentation_pages/selections.html)
#replicate_network_analysis calls this, to generate the selections of water molecules in a bin of interest to use in hydrogen bonding analysis
#############
def create_residue_selection(wat_res_num_list):

    #If there are no water molecules, return None
    if (len(wat_res_num_list) == 0):
        return None

    #If there is one molecule, return the selection with only this water molecule
    if (len(wat_res_num_list) == 1): #Corresponds to a list w/one water molecule
        return f"byres (resname TIP3 and resid {wat_res_num_list[0]})"

    #If there are two or more water molecules, return the selection with all water molecules of interest
    if (len(wat_res_num_list) > 1):

        #Set up selection with byres syntax and resname and first water resid
        water_selection = f"byres (resname TIP3 and (resid {wat_res_num_list[0]}"

        #Add all of the other water molecules besides the first one to the selection
        for r in wat_res_num_list[1:]:
            water_selection = water_selection + f" or resid {r}"

        #Add closing parentheses to the selection and return
        water_selection = water_selection + "))"
        return water_selection

############
#Channel axis projection code written by Michiel Niesen. Many thanks for sharing!
############
def channelaxis(dirv,cv): #thanks Michiel for this function!
    '''
    Determine a normalized vector along the channel axis, based on a direction and a center point.
    '''
    r = np.dot(dirv,cv)
    if r<0:
        return -1*cv/npla.norm(cv)
    else:
        return cv/npla.norm(cv)

############
#Run MDAnalysis call to obtain projection of a specific residue for a specific structure on the channel axis- this will be either for one frame of a trajectory (calc_channel_projection will then include the output in an array for the trajectory) or for a single structure, called by calc_channel_projection function (single structure analysis does not occur in this script)
############
def mdanalysis_individual_projection_value(center, direction, calphas, residue_selection_to_place):
    ############
    #Obtain the channel axis, many thanks to Michiel for this code and above channelaxis function!
    ############
    cc = center.centroid() # Center of the channel
    dc = direction.centroid() - cc # Direction of the channel axis vector
    p1, p2, p3 = calphas.principal_axes() # Principal axes of all Calpha atoms
    caxis = channelaxis(dc,p3) # Calculate normalized vector along channel axis

    ############
    #Obtain the centroid of the four relevant CAs and project along the channel axis
    ############
    residue_centroid = residue_selection_to_place.centroid() - cc # Coordinate of this centroid- 4 CA atoms- relative to the center as defined by Ala17 CA atom centroid, thanks so much Michiel!
    residue_centroid_proj = np.dot(residue_centroid, caxis) # Projection on the channel axis, thanks so much Michiel!

    return residue_centroid_proj

############
#Set up MDAnalysis call to obtain projections of each residue on the channel axis
#Many many thanks to Michiel Niesen for sharing!
#This can be run for a trajectory (is_trajectory argument set to True, from replicate_network_analysis) or a static structure (is_trajectory argument set to False, not used in this script but was used in minimum distance script)
############
def calc_channel_projection(univ_to_analyze_proj, residue_number, args, is_trajectory):

    ############
    #Set up group selections, many thanks to Michiel Niesen for this code
    ############
    center = univ_to_analyze_proj.select_atoms('protein and resname ALA and resid 17 and name CA') # ALA 17 Calphas
    direction = univ_to_analyze_proj.select_atoms('protein and resname TRP and resid 23 and name CA') # TRP 23 Calphas
    calphas = univ_to_analyze_proj.select_atoms('protein and name CA') # All Calphas
    residue_selection_to_place = univ_to_analyze_proj.select_atoms(f"protein and resid {residue_number} and name CA") #Calphas of the residue of interest

    ############
    #If a trajectory is being analyzed, obtain projections for at each time desired (starting at 30 ns) by iterating through the trajectory
    ############
    if is_trajectory:

        time_list = [] #List for times
        proj_list = [] #List for projections

        #Iterate through the trajectory- at each point obtain the location of the relevant residue with respect to the channel axis
        for f in univ_to_analyze_proj.trajectory[3000::args.frame_stride]: #Analyze the trajectory starting at 30 ns (consulted https://github.com/MDAnalysis/MDAnalysisCookbook/blob/master/examples/blocks.py, https://www.mdanalysis.org/MDAnalysisTutorial/trajectories.html?highlight=select%20frames has stride)

            time_list.append(univ_to_analyze_proj.trajectory.time) #Obtain time and add to list (https://www.mdanalysis.org/MDAnalysisTutorial/trajectories.html?highlight=trajectory%20times)
            new_residue_centroid_proj = mdanalysis_individual_projection_value(center, direction, calphas, residue_selection_to_place) #Obtain projection of this residue for this time
            proj_list.append(new_residue_centroid_proj) #Add on the projection to the growing projection list

        #Change lists to arrays after population, for efficiency (https://stackoverflow.com/questions/7332841/add-single-element-to-array-in-numpy)
        time_array = np.array(time_list)
        proj_array = np.array(proj_list)

        return proj_array, time_array #Return the projection and time arrays

    ############
    #If a structure is being analyzed, obtain projection for this structure- because there is one structure, no iteration is carried out
    ############
    if not is_trajectory:
        residue_centroid_proj_value = mdanalysis_individual_projection_value(center, direction, calphas, residue_selection_to_place) #Obtain projection for this residue
        return residue_centroid_proj_value #Return this projection

############
#Carry out the H bond network analysis for each replicate
############
def replicate_network_analysis(args, ch_ax_list, replicate_directory, analysis_output):

    #############
    #Set up- load universe, renumber residues, create data structures
    #############
    #Load the universe for the replicate specified in the function input replicate_directory
    input_xtc = f"{replicate_directory}/md_1.xtc"
    input_tpr = f"{replicate_directory}/md_1.tpr"
    rep_name = replicate_directory.split("/")[-1] #Name for output purposes
    analysis_output.write(f"Values for {rep_name}\n")
    univ_for_analysis = MDAnalysis.Universe(input_tpr, input_xtc)

    #Renumber the residues of the simulation universe so that all four chains have the same residue numbering. Users should check their universes to see whether or not any or all need renumbering.
    residue_numbers = np.array([at.resid for at in univ_for_analysis.select_atoms("protein")]) #Obtain all residue numbers
    residue_numbers_condensed_list = [r for r, i in groupby(residue_numbers)] #Remove duplicates- this gives 1 resid/residue (https://stackoverflow.com/questions/38065898/how-to-remove-the-adjacent-duplicate-value-in-a-numpy-array)
    residue_numbers_condensed = np.array(residue_numbers_condensed_list)
    residue_numbers_condensed_incremented = residue_numbers_condensed + 1 #Increment protein indices by 1
    residue_numbers_condensed_incremented_mod = residue_numbers_condensed_incremented % 33 #Obtain indices mod 33
    residue_numbers_condensed_incremented_mod[residue_numbers_condensed_incremented_mod == 0] = 33 #Replace 0 with 33 (https://www.w3resource.com/python-exercises/numpy/python-numpy-exercise-88.php, https://stackoverflow.com/questions/3803331/how-can-i-modulo-when-my-numbers-start-from-1-not-zero)
    univ_for_analysis.select_atoms("protein").residues.resids = residue_numbers_condensed_incremented_mod #Set residue IDs to the desired values in the simulation universe

    #Set up bottleneck count list, to be populated with each frame's bottleneck count
    bottleneck_counts = [] #Count of bottlenecks at each time

    #Set up dictionaries for hydrogen bonding statistics. Each has keys of all channel axis slices (except the most C-terminal one, since that is only included as being analyzed as adjacent to the next-most-C-terminal one) and values of lists of the relevant metric for each frame analyzed
    hb_dictionary = dict.fromkeys(ch_ax_list[0:-1]) #Hydrogen bond counts between water molecules at each channel axis slice and those at the next-most-C-terminal slice
    hb_diff_dictionary = dict.fromkeys(ch_ax_list[0:-1]) #Hydrogen bond directionality metrics between water molecules at each channel axis slice and those at the next-most-C-terminal slice
    bottlenecks_dictionary = dict.fromkeys(ch_ax_list[0:-1]) #Presence or absence of a bottleneck between water molecules at each channel axis slice and those at the next-most-C-terminal slice (the lists will have either 0 or 1 as each entry)

    #Initialize empty lists
    for ch_ax_k in hb_dictionary.keys():
        hb_dictionary[ch_ax_k] = []
        hb_diff_dictionary[ch_ax_k] = []
        bottlenecks_dictionary[ch_ax_k] = []

    frame_count = 0 #Frame counter- for reporting to see if replicates and simulation conditions have the same frame counts
    cylR2 = 225.0 #Maximum distance, squared, from the channel axis for waters to be included in the analysis

    #Set up selections, many thanks to Michiel Niesen for sharing this code!
    center = univ_for_analysis.select_atoms('protein and resname ALA and resid 17 and name CA') #ALA17 Calphas
    direction = univ_for_analysis.select_atoms('protein and resname TRP and resid 23 and name CA') #TRP23 Calphas
    calphas = univ_for_analysis.select_atoms('protein and name CA') #All Calphas
    waters = univ_for_analysis.select_atoms('resname TIP3 and not (name H* or name [123]H or type H)') #Water oxygens

    #Set up channel binning
    increment_for_edges = float(args.axis_increments) / 2.0 #define increment for edges- to ensure consistent bin size throughout, ** This is an update and was not used in the script to generate the Communications Biology figures. Instead, the channel axis bounds were args.minimum_axis_value and args.maximum_axis_value: so the bins at the edges were half as large as the other bins. With this update, all bins should be the same size.
    minimum_axis_value_for_inclusion = args.minimum_axis_value - increment_for_edges #minimum axis projection for inclusion
    maximum_axis_value_for_inclusion = args.maximum_axis_value + increment_for_edges #maximum axis projection for inclusion
    analysis_output.write(f"Minimum projection for inclusion {minimum_axis_value_for_inclusion}\nMaximum projection for inclusion {maximum_axis_value_for_inclusion}\n")

    #############
    #If the user wants to label residues with their channel axis projections, calculate the channel axis projections and determine the label for each residue
    #Note this involves a sacrifice of runtime for modularity- some calculations in calc_channel_projection are repeated in the H bond analysis code below, but this arrangement enables using the calc_channel_projection code used in another script, and thus this arrangement was selected
    #############
    if (args.residue_numbers_for_label != None):

        #Set up dictionaries and time checking list for H bond analysis that will be run
        ch_proj_dictionary = dict.fromkeys(args.residue_numbers_for_label) #Dictionary with keys of residue numbers whose projections will be labelled and values of np arrays corresponding to projections for this residue at each time point analyzed for this replicate
        name_num_dictionary = dict.fromkeys(args.residue_numbers_for_label) #Dictionary with keys of residue numbers whose projections will be labelled and values of labels of the form (resname)(resnumber)
        times_for_h_bond_calc_list = [] #List of times used in H bond analysis, to compare to times used in projection calculations

        #For each residue, calculate the projections and determine the label
        for k_of_res in ch_proj_dictionary:
            ch_proj_dictionary[k_of_res], times_for_projection_calc = calc_channel_projection(univ_for_analysis, k_of_res, args, True) #Each residue key has a value of a np array that is the projection array. Also record the projection times to check against the H bond analysis times and confirm they are identical. (Analyses for each residue should yield identical times, so times_for_projection_calc will be overwritten, which should not cause problems.)

            #Convert the residue number to a label of the form (resname)(resnumber)
            res_name_for_residue_letter_number_ID = univ_for_analysis.select_atoms(f"protein and resid {k_of_res}")[0].resname #This takes the resname for the first atom in the selection corresponding to this residue: because there is one residue number selected, all atoms in the selection should have the same resname
            residue_letter_number_ID = f"{res_name_for_residue_letter_number_ID}{k_of_res}"
            #If HSD or HSP is in the name, replace with HIS (https://stackoverflow.com/questions/3437059/does-python-have-a-string-contains-substring-method)
            if (("HSD" not in residue_letter_number_ID) and ("HSP" not in residue_letter_number_ID)):
                residue_letter_number_ID_to_use = residue_letter_number_ID
            if ("HSD" in residue_letter_number_ID):
                residue_letter_number_ID_to_use = residue_letter_number_ID.replace("HSD", "HIS")
            if ("HSP" in residue_letter_number_ID):
                residue_letter_number_ID_to_use = residue_letter_number_ID.replace("HSP", "HIS")
            name_num_dictionary[k_of_res] = residue_letter_number_ID_to_use

    #If the user does not want projection labeling, set the dictionaries used to store projections and residue labels to be None
    if (args.residue_numbers_for_label == None):
        ch_proj_dictionary = None
        name_num_dictionary = None
    
    #############
    #Run hydrogen bonding analysis for each trajectory frame of interest
    #############
    for f in univ_for_analysis.trajectory[3000 :: args.frame_stride]: #Analyze the trajectory starting at 30 ns (consulted https://github.com/MDAnalysis/MDAnalysisCookbook/blob/master/examples/blocks.py, https://www.mdanalysis.org/MDAnalysisTutorial/trajectories.html?highlight=select%20frames has stride)

        #############
        #Extract frame/time information
        #############
        frame_to_use = univ_for_analysis.trajectory.frame #Extract frame number
        wat_resid_time = univ_for_analysis.trajectory.time #Extract time
        frame_count += 1 #Increment frame count by 1

        #If the user wants to label residues with their channel axis projections, record the H bonding analysis times in the list times_for_h_bond_calc_list, to compare to times_for_projection_calc and confirm they are equal, after all frames have been iterated over
        if (args.residue_numbers_for_label != None):
            times_for_h_bond_calc_list.append(wat_resid_time)

        #############
        #Set up dictionary and counter to hold information for waters for this frame only
        #############
        #Dictionary with keys of channel axis coordinates, and values of lists of water molecules binned for each channel axis coordinate for this frame
        c_ax_dictionary = dict.fromkeys(ch_ax_list)
        for c_ax in c_ax_dictionary.keys():

            #Initialize values to be empty lists, will be populated with water molecules
            c_ax_dictionary[c_ax] = []

        #Counter to count how many bottlenecks there are for this frame, initialize to 0, will increment when a bottleneck is found
        bottleneck_counter = 0

        ###########
        #For each water, determine its channel axis coordinate and how far it is from the channel axis
        #Many thanks to Michiel Niesen for sharing this binning code!
        ###########
        #Define selections for binning calculation
        cc = center.centroid() # Center of the channel
        dc = direction.centroid() - cc # Direction of the channel axis vector
        p1, p2, p3 = calphas.principal_axes() # Principal axes of all Calpha atoms
        caxis = channelaxis(dc,p3) # Calculate normalized vector along channel axis

        #Check each water- if it is within 15 A of the channel axis and within the desired channel axis coordinate values, find the channel axis coordinate bin to which it is closest and add the water to the list for that bin in the c_ax_dictionary dictionary list entry
        for wi in range(len(waters)):
            posO = waters[wi].position - cc #Distance vector between oxygen and channel center
            posOx = np.dot(posO,caxis) #Projection of oxygen onto channel coordinate axis
            posOyz2 = npla.norm(posO)**2 - abs(posOx)**2 #This calculates how far the oxygen is from the channel axis coordinate line

            #If the water is within the desired channel axis coordinate bounds, and it is within square root of cylR2 from the channel axis coordinate line, determine the appropriate bin and place this water in the list corresponding to that bin in c_ax_dictionary
            if ((posOx > minimum_axis_value_for_inclusion) and (posOx < maximum_axis_value_for_inclusion) and (posOyz2 < cylR2)):
                relevant_c_ax_coord = min(ch_ax_list, key = lambda x : abs(x - posOx)) #Determine to which coordinate bin this water is the closest (code line from https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value, this checks the different channel axis coordinates used in binning to find the one with the smallest difference between itself and posOx, giving the channel axis coordinate value used for binning closest to posOx)
                c_ax_dictionary[relevant_c_ax_coord].append(waters[wi]) #Add the water to the appropriate channel axis coordinate bin

        ###########
        #Iterate through all channel axis coordinate bins, running the hydrogen bonding analyses to examine hydrogen bonds between waters in the bin being examined and that in the immediately adjacent bin closer to the C-terminus (because the next most C-terminal bin is also considered in each iteration, the last/most C-terminal entry in ch_ax_list is not going to have an iteration for itself in this for loop)
        ###########
        for coord in ch_ax_list[:-1]:             
            #Obtain resids of waters in this channel axis bin (first_resid_list) and those in the immediately adjacent bin closer to the C-terminus (second_resid_list)
            #Caution- if the axis increments/minimum channel axis coordinate/maximum channel axis coordinate are not carefully specified, adding args.axis_increments to each bin coordinate in the dictionary to obtain the adjacent one may not work
            first_resid_list = [water_molecule.resid for water_molecule in c_ax_dictionary[coord]]
            second_resid_list = [water_molecule.resid for water_molecule in c_ax_dictionary[coord + args.axis_increments]]

            #Set up selections for the MDAnalysis.analysis.hbonds.HydrogenBondAnalysis call
            first_resid_selection = create_residue_selection(first_resid_list) #Set up selections
            second_resid_selection = create_residue_selection(second_resid_list)

            ###########
            #Run H bond analysis, considering H bonds between water molecules in this channel axis coordinate bin and those in the immediately adjacent bin closer to the C-terminus (https://www.mdanalysis.org/docs/documentation_pages/analysis/hbond_analysis.html)
            #Run two MDAnalysis.analysis.hbonds.HydrogenBondAnalysis calls: one in which the first selection is the acceptor, and one in which the first selection is the donor
            ###########
            if (first_resid_selection != None and second_resid_selection != None): #Only run the analysis if there is at least one atom in both selections

                #Determine count of H bonds between water molecules in this channel axis coordinate bin and those in the immediately adjacent bin closer to the C-terminus in which those waters closer to the N-terminus are acceptors/those closer to the C-terminus are donors
                hbanalysis_acceptors = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(univ_for_analysis, first_resid_selection, second_resid_selection, selection1_type = "acceptor")
                hbanalysis_acceptors.run(start = frame_to_use, stop = frame_to_use + 1, selection1_type = "acceptor") #Only run the analysis for the frame currently being analyzed
                h_b_time_a = hbanalysis_acceptors.count_by_time()[0][0] #Find the time to confirm the correct frame was in fact used- uses the fact there is only one frame in the analysis (details on the count_by_time() output are found here: https://www.mdanalysis.org/docs/documentation_pages/analysis/hbond_analysis.html)
                h_b_count_a = hbanalysis_acceptors.count_by_time()[0][1] #Find how many hydrogen bonds there are- uses the fact there is only one frame in the analysis (details on the count_by_time() output are found here: https://www.mdanalysis.org/docs/documentation_pages/analysis/hbond_analysis.html)

                #Determine count of H bonds between water molecules in this channel axis coordinate bin and those in the immediately adjacent bin closer to the C-terminus in which those waters closer to the N-terminus are donors/those closer to the C-terminus are acceptors
                hbanalysis_donors = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(univ_for_analysis, first_resid_selection, second_resid_selection, selection1_type = "donor")
                hbanalysis_donors.run(start = frame_to_use, stop = frame_to_use + 1, selection1_type = "donor") #Only run the analysis for the frame currently being analyzed
                h_b_time_d = hbanalysis_donors.count_by_time()[0][0] #Find the time to confirm the correct frame was in fact used- uses the fact there is only one frame in the analysis (details on the count_by_time() output are found here: https://www.mdanalysis.org/docs/documentation_pages/analysis/hbond_analysis.html)
                h_b_count_d = hbanalysis_donors.count_by_time()[0][1] #Find how many hydrogen bonds there are- uses the fact there is only one frame in the analysis (details on the count_by_time() output are found here: https://www.mdanalysis.org/docs/documentation_pages/analysis/hbond_analysis.html)

                #Sum donor and acceptor H bond counts to obtain total H bond count
                h_b_count = h_b_count_a + h_b_count_d

                #Extra-cautious time safety checks to confirm correct frame was analyzed
                if (h_b_time_a != wat_resid_time):
                    warnings.warn("***Times indicate different frames pulled for water ID and H bond (selection1 being acceptors) calculations")
                    analysis_output.write("***Times indicate different frames pulled for water ID and H bond (selection1 being acceptors) calculations\n")

                if (h_b_time_d != wat_resid_time):
                    warnings.warn("***Times indicate different frames pulled for water ID and H bond (selection1 being donors) calculations") 
                    analysis_output.write("***Times indicate different frames pulled for water ID and H bond (selection1 being donors) calculations\n")

            #If there are no waters in one or both selections, then there are no hydrogen bonds between water molecules at the relevant channel axis slices, so all hydrogen bond counts are 0
            else:
                h_b_count = 0
                h_b_count_a = 0
                h_b_count_d = 0

            #Update dictionaries with data for this channel axis coordinate being analyzed
            #Add the total H bond count to the list for the relevant channel axis coordinate slice in hb_dictionary
            hb_dictionary[coord].append(h_b_count)

            #Add the H bond directionality metric to the list for the relevant channel axis coordinate slice in hb_diff_dictionary (count of H bonds between water molecules in this channel axis coordinate bin and those in the immediately adjacent bin closer to the C-terminus in which those waters closer to the N-terminus are acceptors minus count of H bonds between water molecules in this channel axis coordinate bin and those in the immediately adjacent bin closer to the C-terminus in which those waters closer to the N-terminus are donors)
            hb_diff_dictionary[coord].append(h_b_count_a - h_b_count_d)

            #If there is a bottleneck, increment the bottleneck counter and record the presence of a bottleneck at the relevant channel axis coordinate slice in bottlenecks_dictionary (by appending a 1 to the appropriate list). If there is not a bottleneck, record the absence of a bottleneck at the relevant channel axis coordinate slice in bottlenecks_dictionary (by appending a 0 to the appropriate list).
            if (h_b_count == 0):
                bottlenecks_dictionary[coord].append(1)
                bottleneck_counter += 1 #If there are no H bonds between water molecules in this channel axis coordinate bin and those in the immediately adjacent bin closer to the C-terminus, increment the bottleneck counter to record the existence of a bottleneck for this slice for this frame
            else:
                bottlenecks_dictionary[coord].append(0)
            
        #Add count of bottlenecks for this frame to the list of frame-level bottleneck counts (this is done once for each frame, after all MDAnalysis calls' results have beeen recorded)
        bottleneck_counts.append(bottleneck_counter)

    #If projection calculations were run as well, times for these calculations and the H bonding calculations should be equal- this is a check to confirm the H bonding metric data and the projection data correspond to the same frames
    if (args.residue_numbers_for_label != None):
        if (not np.array_equal(np.array(times_for_h_bond_calc_list), times_for_projection_calc)): #Use np to confirm the times are all equal (https://docs.scipy.org/doc/numpy/reference/generated/numpy.array_equal.html)
            warnings.warn("***Times indicate different frames pulled for hydrogen bonding and projection calculations")
            analysis_output.write("***Times indicate different frames pulled for hydrogen bonding and projection calculations\n")

    #Change lists to arrays after population, for efficiency (https://stackoverflow.com/questions/7332841/add-single-element-to-array-in-numpy)
    bottleneck_counts_array = np.array(bottleneck_counts)
    for ch_ax_key in hb_dictionary.keys():
        hb_dictionary[ch_ax_key] = np.array(hb_dictionary[ch_ax_key])
        hb_diff_dictionary[ch_ax_key] = np.array(hb_diff_dictionary[ch_ax_key])
        bottlenecks_dictionary[ch_ax_key] = np.array(bottlenecks_dictionary[ch_ax_key])

    #Report bottlenecks across channel
    analysis_output.write(f"Average bottleneck count across channel : {np.average(bottleneck_counts_array)}, obs count: {frame_count}, sem: {scipy.stats.sem(bottleneck_counts_array)}\n")

    #Return np array of bottleneck counts; dictionaries of H bond count, H bond directionality metric, and bottleneck presence/absence for each channel axis slice; count of frames analyzed; dictionary of channel axis projections of residues if wanted for plotting; and residue number/label mapping dictionary if wanted for plotting
    return bottleneck_counts_array, hb_dictionary, hb_diff_dictionary, bottlenecks_dictionary, frame_count, ch_proj_dictionary, name_num_dictionary

############
#Generate histogram of bottleneck counts
############
def bottleneck_histogram_plot(pooled_dictionary_of_analyses, args, analysis_output, dictionary_of_colors_for_plotting):

    ############
    #Generate histogram
    ############
    #Determine the maximum bottleneck count to determine the largest bin value (iterate through the dictionary and update running_max [initially set to 0] each time a bottleneck count larger than any seen previously is seen)
    running_max = 0
    for kb in pooled_dictionary_of_analyses.keys():
        max_bottlenecks_for_key = max(pooled_dictionary_of_analyses[kb]["bottleneck_count_array"])
        if (running_max < max_bottlenecks_for_key):
            running_max = max_bottlenecks_for_key

    #Use this maximum bottleneck count to generate bins for the histogram (https://matplotlib.org/3.2.1/api/_as_gen/matplotlib.pyplot.hist.html)
    bins_for_hist = np.arange(0, running_max + 3, 1)

    #Plot histogram for each simulation condition, also record bin populations in output file (https://matplotlib.org/3.2.1/api/_as_gen/matplotlib.pyplot.hist.html)
    for sim_key in pooled_dictionary_of_analyses.keys(): #Iterate over simulation conditions
        analysis_output.write("HISTOGRAM\n")
        analysis_output.write(f"{sim_key}\n")
        val, ed, pat = plt.hist(pooled_dictionary_of_analyses[sim_key]["bottleneck_count_array"], bins = bins_for_hist, align = "left", color = dictionary_of_colors_for_plotting[sim_key]["data"], label = sim_key, alpha = 0.50)
        for bin_val, count_at_bin in zip(ed, val):
            analysis_output.write(f"{bin_val} bottleneck count seen in {count_at_bin} frames, {count_at_bin * 100 / len(pooled_dictionary_of_analyses[sim_key]['bottleneck_count_array'])}% of all frames\n")

    ############
    #Label/format axes, set axis title/label font size, assign the plot a title, and save it (https://stackoverflow.com/questions/9651092/my-matplotlib-pyplot-legend-is-being-cut-off/42303455)
    ############
    #Add title, label axes, and set font size
    plt.title(f"Bottleneck Count\nJob {args.title_of_job}", fontsize = 18)
    plt.xlabel("Bottleneck Count", fontsize = 24)
    plt.ylabel("Frequency", fontsize = 20)
    plt.xticks(fontsize = 16) #Change axis label size (https://www.science-emergence.com/Articles/How-to-change-the-size-of-axis-labels-in-matplotlib-/)
    plt.yticks(fontsize = 22) #Change axis label size (https://www.science-emergence.com/Articles/How-to-change-the-size-of-axis-labels-in-matplotlib-/)

    #Add legend, consulted https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot/27355247 for legend help
    plt.legend(loc = "upper left", bbox_to_anchor = (1.02, 1.00)) #Consulted https://matplotlib.org/tutorials/intermediate/legend_guide.html
    plt.tight_layout() #Set layout (https://matplotlib.org/3.2.1/tutorials/intermediate/tight_layout_guide.html)

    #Save figure
    plt.savefig(f"{args.directory_for_output}/Bottleneck_Histogram_{args.title_of_job}.png")
    plt.savefig(f"{args.directory_for_output}/Bottleneck_Histogram_{args.title_of_job}.pdf")
    plt.savefig(f"{args.directory_for_output}/Bottleneck_Histogram_{args.title_of_job}.eps") #Note a message "The PostScript backend does not support transparency; partially transparent artists will be rendered opaque." could appear, workarounds can be found here: https://stackoverflow.com/questions/19638773/matplotlib-plots-lose-transparency-when-saving-as-ps-eps

    #Clear/close in preparation for next plot (https://stackoverflow.com/questions/17106288/matplotlib-pyplot-will-not-forget-previous-plots-how-can-i-flush-refresh)
    plt.clf()
    plt.cla()
    plt.close()

############
#Generate line plots of channel axis coordinate versus H bond count, H bond directionality metric, and bottleneck percentage
############
def hb_line_plot(pooled_dictionary_of_analyses, selection_for_dictionary, desc_for_plot, args, analysis_output, dictionary_of_colors_for_plotting, ch_ax_coord_list, res_proj_for_plot, num_name_linking_dictionary):

    #Obtain key names- note this is currently hard-coded to work with two simulation conditions. Edits would be needed to generalize to more simulation conditions.
    first_k = sorted(pooled_dictionary_of_analyses.keys())[0]
    second_k = sorted(pooled_dictionary_of_analyses.keys())[1]

    ############
    #Generate the line plot w/error bars
    ############

    #Obtain sorted list of channel axis coordinate slices (should be identical for both, though extracting separately here)
    first_k_channel_sorted = sorted(pooled_dictionary_of_analyses[first_k][selection_for_dictionary].keys())
    second_k_channel_sorted = sorted(pooled_dictionary_of_analyses[second_k][selection_for_dictionary].keys())

    #Obtain plot averages
    #[np.average(pooled_dictionary_of_analyses[first_k][selection_for_dictionary][ch_a_1]) for ch_a_1 in first_k_channel_sorted] will give the average for each channel axis coordinate as a list in the desired order (sorting is needed here and in SEM calculation for order)
    first_k_averages = [np.average(pooled_dictionary_of_analyses[first_k][selection_for_dictionary][ch_a_1]) for ch_a_1 in first_k_channel_sorted]
    second_k_averages = [np.average(pooled_dictionary_of_analyses[second_k][selection_for_dictionary][ch_a_2]) for ch_a_2 in second_k_channel_sorted]

    #Determine the error bars- for H bond count and directionality metric, obtain for the first_k condition and second_k condition the SEM for the replicate-level averaged value for each channel axis coordinate slice
    #[scipy.stats.sem(pooled_dictionary_of_analyses[first_k][selection_for_dictionary][ch_a_1]) for ch_a_1 in first_k_channel_sorted] will give the SEM for each channel axis coordinate as a list in the desired order
    if (selection_for_dictionary != "overall_bottleneck_dictionary_by_ch_ax"):
        first_k_error = [scipy.stats.sem(pooled_dictionary_of_analyses[first_k][selection_for_dictionary][ch_a_1]) for ch_a_1 in first_k_channel_sorted]
        second_k_error = [scipy.stats.sem(pooled_dictionary_of_analyses[second_k][selection_for_dictionary][ch_a_2]) for ch_a_2 in second_k_channel_sorted]

    #Determine the error bars-for bottleneck percentage, the error bar is 0, since this is a percentage of frames (https://numpy.org/doc/stable/reference/generated/numpy.zeros.html)
    else:
        first_k_error = np.zeros(len(first_k_channel_sorted))
        second_k_error = np.zeros(len(second_k_channel_sorted))

    #For each simulation condition, plot the averaged value for each channel axis coordinate, along with error bars
    plt.errorbar(first_k_averages, first_k_channel_sorted, color = dictionary_of_colors_for_plotting[first_k]["data"], xerr = first_k_error, label = first_k, ecolor = dictionary_of_colors_for_plotting[first_k]["uncertainty"], markersize = 5)
    plt.errorbar(second_k_averages, second_k_channel_sorted, color = dictionary_of_colors_for_plotting[second_k]["data"], xerr = second_k_error, label = second_k, ecolor = dictionary_of_colors_for_plotting[second_k]["uncertainty"], markersize = 5)

    #Record the data plotting
    analysis_output.write("LINE PLOT\n")
    analysis_output.write(f"{selection_for_dictionary} for {first_k}\n")
    for xv_1, yv_1, errv_1 in zip(first_k_averages, first_k_channel_sorted, first_k_error):
        analysis_output.write(f"x {xv_1}, y {yv_1}, errv {errv_1}\n")
    analysis_output.write(f"{selection_for_dictionary} for {second_k}\n")
    for xv_2, yv_2, errv_2 in zip(second_k_averages, second_k_channel_sorted, second_k_error):
        analysis_output.write(f"x {xv_2}, y {yv_2}, errv {errv_2}\n")
    analysis_output.write("\n")

    ############
    #If user wants to label certain residues- obtain the projected average coordinate value, draw a line, label
    ############
    #Set axis limits if wanted (https://matplotlib.org/3.2.2/api/_as_gen/matplotlib.pyplot.ylim.html?highlight=ylim#matplotlib.pyplot.ylim)
    if (args.minimum_y_for_plot != None):
        plt.ylim(bottom = args.minimum_y_for_plot)
    if (args.maximum_y_for_plot != None):
        plt.ylim(top = args.maximum_y_for_plot)

    #If wanted, add on dashed lines for the relevant residues, draw a line for each residue of interest corresponding to its average projection
    if (args.residue_numbers_for_label != None):
        minimum_x_value, maximum_x_value = plt.xlim() #Obtain axis limits (https://matplotlib.org/3.2.2/api/_as_gen/matplotlib.pyplot.xlim.html?highlight=xlim#matplotlib.pyplot.xlim)
        close_to_max = minimum_x_value + 0.85 * (maximum_x_value - minimum_x_value) #Obtain x coordinate near maximum

        #Label projection of each residue
        for res_num_to_mark in args.residue_numbers_for_label:
            projected_coord_of_res = res_proj_for_plot[res_num_to_mark] #Obtain projected-onto-channel-axis value for residue of interest from the dictionary

            #Draw a horizontal line at the channel axis coordinate corresponding to this projection (https://matplotlib.org/3.2.2/api/_as_gen/matplotlib.pyplot.axhline.html?highlight=axhline#matplotlib.pyplot.axhline), and have the line width be small (https:/matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D.set_linewidth)
            plt.axhline(projected_coord_of_res, linestyle = "--", linewidth = 1, color = (0.07, 0.61, 1.00))

            #Write the residue number and associated projected channel axis coordinate value to the output file
            analysis_output.write(f"PROJECTION for {res_num_to_mark} {projected_coord_of_res}\n")

            #Add text to associate the appropriate residue name/number identifier with the projection line just drawn (offsets added to x axis minimum and projected_coord_of_res for cleanliness), using num_name_linking_dictionary to find the correct residue name/number identifier
            plt.text(close_to_max, (projected_coord_of_res - 0.25), num_name_linking_dictionary[res_num_to_mark], size = 11, ha = "left", color = (0.07, 0.61, 1.00)) #Add text (https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/arrow_demo.html, https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.text.html)

    #For H bond count plot generate a vertical line at 0 to facilitate interpretation, many thanks to Martin Gelenter and Venkata Shiva Mandala for this suggestion! (https://matplotlib.org/3.2.2/api/_as_gen/matplotlib.pyplot.axvline.html?highlight=axvline#matplotlib.pyplot.axvline)
    if (selection_for_dictionary == "hb_dictionary_by_ch_ax"):
        plt.axvline(0.0, linestyle = "--", linewidth = 1, color = (0.47, 0.29, 0.00)) 

    ############
    #Label/format axes, set axis title/label font size, assign the plot a title, and save it (https://stackoverflow.com/questions/9651092/my-matplotlib-pyplot-legend-is-being-cut-off/42303455)
    ############
    #Add title, label axes and set font size, also flip y axis
    plt.title(f"{desc_for_plot} Line Plot\nJob {args.title_of_job}", fontsize = 18)
    plt.ylabel("Channel Axis Coordinate($\AA$)", fontsize = 24) #Consulted https://idlangstrom.wordpress.com/2014/12/05/angstrom-and-other-astronomical-symbols-in-python/ for angstrom formatting
    plt.xlabel(desc_for_plot, fontsize = 20)
    plt.xticks(fontsize = 14) #Change axis label size (https://www.science-emergence.com/Articles/How-to-change-the-size-of-axis-labels-in-matplotlib-/)
    plt.yticks(fontsize = 22) #Change axis label size (https://www.science-emergence.com/Articles/How-to-change-the-size-of-axis-labels-in-matplotlib-/)
    plt.gca().invert_yaxis() #Flip y axis (https://stackoverflow.com/questions/2051744/reverse-y-axis-in-pyplot)

    #Add legend, consulted https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot/27355247 for legend help
    plt.legend(loc = "upper left", bbox_to_anchor = (1.02, 1.00)) #Consulted https://matplotlib.org/tutorials/intermediate/legend_guide.html
    plt.tight_layout() #Set layout (https://matplotlib.org/3.2.1/tutorials/intermediate/tight_layout_guide.html)

    #Save figure
    plt.savefig(f"{args.directory_for_output}/{selection_for_dictionary}_Line_Chart_{args.title_of_job}.png")
    plt.savefig(f"{args.directory_for_output}/{selection_for_dictionary}_Line_Chart_{args.title_of_job}.pdf")
    plt.savefig(f"{args.directory_for_output}/{selection_for_dictionary}_Line_Chart_{args.title_of_job}.eps") #Note a message "The PostScript backend does not support transparency; partially transparent artists will be rendered opaque." could appear, workarounds can be found here: https://stackoverflow.com/questions/19638773/matplotlib-plots-lose-transparency-when-saving-as-ps-eps

    #Clear/close in preparation for next plot (https://stackoverflow.com/questions/17106288/matplotlib-pyplot-will-not-forget-previous-plots-how-can-i-flush-refresh)
    plt.clf()
    plt.cla()
    plt.close()

if __name__ == "__main__":

    #############
    #Parse arguments, make output directory/files, load relevant dictionaries, set up channel axis bins
    #############
    #Parse arguments
    parser = argparse.ArgumentParser(description = "Plots metrics for hydrogen bonds between water molecules in different channel axis coordinate bins")
    parser.add_argument("--directory_for_output", required = True, type = str, help = "Directory to be created to hold the results from this analysis")
    parser.add_argument("--title_of_job", required = True, type = str, help = "Name of job. Best to make distinct for each run.")
    parser.add_argument("--dictionary_of_simulations", required = True, type = str, help = "Path of a file (json format) containing a dictionary for accessing simulation data. Keys are shorthand for simulations. Values are dictionaries themselves- named key of 'parent_directory' corresponds to value of directory stem with .tpr and .xtc files.")
    parser.add_argument("--selected_simulations", required = True, nargs = "+", type = str, help = "Keys in --dictionary_of_simulations json to be used")
    parser.add_argument("--dictionary_of_colors", required = True, type = str, help = "Path of a file (json format) containing a dictionary of colors for plotting. Keys are those in --dictionary_of_simulations. Values are dictionaries themselves- named keys of 'data' and 'uncertainty' correspond each to a list of three floating point numbers representing the colors to use in plotting the values for this set of simulations and the error bars for this set of simulations, respectively.")
    parser.add_argument("--minimum_axis_value", required = True, type = float, help = "Minimum channel axis value for center of a bin binning water molecules (actual bound is also determined by increments, which is updated from the Communications Biology figure)")
    parser.add_argument("--maximum_axis_value", required = True, type = float, help = "Maximum channel axis value for center of a bin binning water molecules (actual bound is also determined by increments, which is updated from the Communications Biology figure)")
    parser.add_argument("--axis_increments", required = True, type = float, help = "Axis increment values for binning water molecules")
    parser.add_argument("--frame_stride", required = True, type = int, help = "Frame stride. This should be small enough that at least two frames are included for each trajectory.")
    parser.add_argument("--minimum_y_for_plot", type = float, help = "Minimum y value for the plot")
    parser.add_argument("--maximum_y_for_plot", type = float, help = "Maximum y value for the plot")
    parser.add_argument("--residue_numbers_for_label", nargs = "+", type = int, help = "Residue numbers for labelling on plot, if wanted")
    args = parser.parse_args()

    #Create output directory and results file (directory cannot already exist: https://www.geeksforgeeks.org/create-a-directory-in-python/)
    os.mkdir(args.directory_for_output)
    analysis_output = open(f"{args.directory_for_output}/Analysis_Output_{args.title_of_job}.txt", "w")

    #Load color dictionary 
    with open(args.dictionary_of_colors) as color_dictionary_file_name:
        dictionary_of_colors_for_plotting = json.load(color_dictionary_file_name)

    #Load the directories for analysis dictionary
    with open(args.dictionary_of_simulations) as information_dictionary_file_name:
        directories_for_analysis = json.load(information_dictionary_file_name)

    #Obtain projected channel axis bins, which will be used to bin water molecules (https://numpy.org/doc/stable/reference/generated/numpy.arange.html)
    ch_ax_list = np.arange(args.minimum_axis_value, args.maximum_axis_value + args.axis_increments, args.axis_increments)

    #############
    #Run water network analysis calculations
    #############

    #Set up dictionary to hold results pooled for each simulation condition's different replicates, with keys of the simulation condition (corresponding to relevant keys in directories_for_analysis) and values of dictionaries of data with named keys storing data from analyses that will be used in plotting (see information on each key below)
    pooled_dictionary_of_analyses = {}

    #Also set up list to confirm same frame counts are being used for each condition in determining average projections
    frame_count_list = []

    #For each set of simulations, run hydrogen bond network analysis
    for selected_simulation_condition in args.selected_simulations: #Analyze each set of simulations
        pooled_dictionary_of_analyses[selected_simulation_condition] = {}
        pooled_dictionary_of_analyses[selected_simulation_condition]["bottleneck_count_array"] = None #Array of bottleneck count for each frame analyzed
        pooled_dictionary_of_analyses[selected_simulation_condition]["total_count_fr"] = None #Total count of frames analyzed across replicates for this simulation condition
        pooled_dictionary_of_analyses[selected_simulation_condition]["num_to_name"] = None #Dictionary with keys of residue numbers to label and values of labels of form (resname)(resnumber) (To be used for labelling the plot if projections are wanted)
        pooled_dictionary_of_analyses[selected_simulation_condition]["residue_avg_proj"] = None #Dictionary with keys of residue numbers to label and values of np array of average projection onto the channel axis coordinate for each replicate (To be used for labelling the plot if projections are wanted)
        pooled_dictionary_of_analyses[selected_simulation_condition]["AV_hb_dictionary_by_ch_ax"] = None #Dictionary with keys of channel axis coordinates and values of np array of average H bond counts between waters in this bin and the adjacent bin in the positive z direction for each replicate
        pooled_dictionary_of_analyses[selected_simulation_condition]["AV_hb_diff_dictionary_by_ch_ax"] = None #Dictionary with keys of channel axis coordinates and values of np array of average H bond directionality metrics between waters in this bin and the adjacent bin in the positive z direction for each replicate
        pooled_dictionary_of_analyses[selected_simulation_condition]["AV_overall_bottleneck_dictionary_by_ch_ax"] = None #Dictionary with keys of channel axis coordinates and values of np array of percentages of frames with bottlenecks between waters in this bin and the adjacent bin in the positive z direction for each replicate 
        pooled_dictionary_of_analyses[selected_simulation_condition]["bottleneck_count_array"], pooled_dictionary_of_analyses[selected_simulation_condition]["total_count_fr"], pooled_dictionary_of_analyses[selected_simulation_condition]["num_to_name"], pooled_dictionary_of_analyses[selected_simulation_condition]["residue_avg_proj"], pooled_dictionary_of_analyses[selected_simulation_condition]["AV_hb_dictionary_by_ch_ax"], pooled_dictionary_of_analyses[selected_simulation_condition]["AV_hb_diff_dictionary_by_ch_ax"], pooled_dictionary_of_analyses[selected_simulation_condition]["AV_overall_bottleneck_dictionary_by_ch_ax"] = network_analysis(selected_simulation_condition, directories_for_analysis, args, ch_ax_list, analysis_output)
        num_dictionary_to_use = pooled_dictionary_of_analyses[selected_simulation_condition]["num_to_name"] #This refreshes for each iteration- because all residues should have the same numbering/naming, using the dictionary obtained in the most recent iteration of the for loop should not be an issue
        frame_count_list.append(pooled_dictionary_of_analyses[selected_simulation_condition]["total_count_fr"])

    #After all simulation conditions have been analyzed, check that all of them have the same count of observations. If they do not, print out a warning that the average projections across conditions will differ from the average across all pooled frames
    number_of_distinct_obs_counts = len(np.unique(frame_count_list)) #Obtain count of unique values (https://kite.com/python/answers/how-to-count-frequency-of-unique-values-in-a-numpy-array-in-python)
    if (number_of_distinct_obs_counts > 1):
         warnings.warn("***Please note observation counts are not the same for all simulation conditions. Be aware that the cross-condition residue projection average, which is reported, will be different from the average obtained by weighting each frame equally.")
         analysis_output.write("***Please note observation counts are not the same for all simulation conditions. Be aware that the cross-condition residue projection average, which is reported, will be different from the average obtained by weighting each frame equally.\n")

    analysis_output.write("---Cross-replicate Projection Averages---\n")

    #Obtain projections of residues of interest onto the channel axis coordinate, averaging across simulation conditions' replicates
    #Please note: condition-level averages, not replicate-level averages, are averaged: so each condition is weighted equally- if different conditions have different numbers of replicates, all replicates will not be weighted equally.
    if (args.residue_numbers_for_label != None):
        overall_proj_of_res = dict.fromkeys(args.residue_numbers_for_label) #Create dictionary with keys of wanted residue numbers, values will be the projections
        for res_to_get_proj in overall_proj_of_res.keys(): #Iterate over the residue numbers- add the projection for each simulation condition to proj_array, and then obtain the average and store in overall_proj_of_res with a key of the residue number and a value of the average projection
            proj_array = np.empty(0) #Empty np array- will hold the projections for each simulation condition
            for simulation_to_analyze in pooled_dictionary_of_analyses.keys(): #For each simulation condition add averaged projection of this residue for different replicates to the array
                proj_array = np.append(proj_array, np.average(pooled_dictionary_of_analyses[simulation_to_analyze]["residue_avg_proj"][res_to_get_proj]))
            overall_proj_of_res[res_to_get_proj] = np.average(proj_array) #In the projections dictionary add the key-value pair of the residue number key and value of average projection across simulation conditions
            analysis_output.write(f"{res_to_get_proj}: Average projection {overall_proj_of_res[res_to_get_proj]} A, SEM {scipy.stats.sem(proj_array)} A\n")

    #If projections are not wanted, set the dictionary to be None
    if (args.residue_numbers_for_label == None):
        overall_proj_of_res = None

    #############
    #Record results and plot
    ############# 
    #Save data (https://stackoverflow.com/questions/30811918/saving-dictionary-of-numpy-arrays). Note from ipython tests when np.load() is called to re-obtain data allow_pickle = True may need to be set. Also note this webpage discusses the format of the load function output.
    data_output_name = f"{args.directory_for_output}/Pooled_dictionary_{args.title_of_job}.npy"
    np.save(data_output_name, pooled_dictionary_of_analyses) 

    #Obtain histogram of bottleneck counts
    bottleneck_histogram_plot(pooled_dictionary_of_analyses, args, analysis_output, dictionary_of_colors_for_plotting)

    #Obtain line plots of the hydrogen bond counts for each channel axis slice, the directionality metric for each channel axis slice, and the bottleneck percentage for each channel axis slice
    hb_line_plot(pooled_dictionary_of_analyses, "AV_hb_dictionary_by_ch_ax", "CR_Average_Num_HB", args, analysis_output, dictionary_of_colors_for_plotting, ch_ax_list, overall_proj_of_res, num_dictionary_to_use) #Hydrogen bond count plot
    hb_line_plot(pooled_dictionary_of_analyses, "AV_hb_diff_dictionary_by_ch_ax", "CR_Average_HB_A_D_Diff", args, analysis_output, dictionary_of_colors_for_plotting, ch_ax_list, overall_proj_of_res, num_dictionary_to_use) #Hydrogen bond directionality metric plot
    hb_line_plot(pooled_dictionary_of_analyses, "AV_overall_bottleneck_dictionary_by_ch_ax", "CR_Percent_Bottlenecks_Across_Fr", args, analysis_output, dictionary_of_colors_for_plotting, ch_ax_list, overall_proj_of_res, num_dictionary_to_use) #Bottleneck percentage plot

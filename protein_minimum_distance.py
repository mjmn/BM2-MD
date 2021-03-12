import argparse
import glob
from itertools import groupby
import json
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis
import numpy as np
import numpy.linalg as npla
import os
import scipy
from scipy import stats
import warnings


#############
#For each replicate of simulation trajectories specified by simulation_selection_name, calculate minimum distances for each residue in which the user is interested
#############
def minimum_distance_for_selection(simulation_selection_name, args, analysis_output, information_dictionary, res_num_for_plot):
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
    #Create dictionaries to hold residue projections and minimum distances for each replicate
    #CR_min_across_reps and CR_proj_across_reps have averages for each replicate. There is no judgment call needed for average/uncertainty calculations across replicates if all trajectories are the same length, as was the case for all trajectories analyzed. If there are trajectories that are of different lengths, then further calculations may want to employ weighting.
    #############
    #To hold averaged minimum distance for each residue for each replicate
    CR_min_across_reps = dict.fromkeys(res_num_for_plot, np.empty(0)) #https://www.programiz.com/python-programming/methods/dictionary/fromkeys

    #To hold averaged projection for each residue for each replicate
    CR_proj_across_reps = dict.fromkeys(res_num_for_plot, np.empty(0)) #https://www.programiz.com/python-programming/methods/dictionary/fromkeys

    #This will hold lengths of each replicate's data array, so a warning can be printed if they differ. Also record total count of frames analyzed.
    CR_lengths = np.empty(0)
    frame_count_total = 0

    #############
    #Run minimum distance and projection analyses for each trajectory, and record the replicate-wide averages in CR_min_across_reps and CR_proj_across_reps, for plotting
    #############
    for replicate_directory in replicate_directories:
        print(f"On {replicate_directory}") #Track progress
        if os.path.isfile(f"{replicate_directory}/md_1.xtc"): #Confirm path exists (https://docs.python.org/dev/library/os.path.html#os.path.isfile)

            #The replicate_distance_and_projection_calc function calculates the minimum distances and projections for each simulation. num_name_map is a dictionary mapping residue numbers to residue number/name codes, a mapping which is the same for each replicate, so its getting overwritten is not an issue. individ_min_rep and individ_proj_rep are dictionaries of minimum distance time series and projection time series (respectively) for each residue for this replicate only, whose averages will be added to CR_min_across_reps and CR_proj_across_reps, respectively.
            num_name_map, individ_min_rep, individ_proj_rep = replicate_distance_and_projection_calc(simulation_selection_name, replicate_directory, args, analysis_output, res_num_for_plot)

            #Record the observation count for this replicate, making use of the fact that all of the residues' distance arrays for a given simulation have the same length so the first residue is used. Also update the frame count total.
            frame_count = len(individ_min_rep[res_num_for_plot[0]])
            CR_lengths = np.append(CR_lengths, frame_count)
            frame_count_total = frame_count_total + frame_count

            #Obtain average minimum distances and projections for each residue for this replicate only, and record in CR_min_across_reps and CR_proj_across_reps
            for res_number_to_add in sorted(res_num_for_plot):
                CR_min_across_reps[res_number_to_add] = np.append(CR_min_across_reps[res_number_to_add], np.average(individ_min_rep[res_number_to_add]))
                CR_proj_across_reps[res_number_to_add] = np.append(CR_proj_across_reps[res_number_to_add], np.average(individ_proj_rep[res_number_to_add]))

    #After all replicates have been analyzed, check that all of them have the same count of observations. If they do not, print out a warning that the average across frames will differ from the average across replicates
    number_of_distinct_obs_counts = len(np.unique(CR_lengths)) #Obtain count of unique values (https://kite.com/python/answers/how-to-count-frequency-of-unique-values-in-a-numpy-array-in-python)
    if (number_of_distinct_obs_counts > 1):
         warnings.warn("***Please note observation counts are not the same for all replicates. Be aware that the cross-replicate average, which is reported, will be different from the average obtained by weighting each frame equally.")
         analysis_output.write("***Please note observation counts are not the same for all replicates. Be aware that the cross-replicate average, which is reported, will be different from the average obtained by weighting each frame equally.\n")

    ############
    #Record averages and SEMs for minimum distances for each residue
    ############
    CR_avg_res_distance_array = np.empty(0) #Array to hold averages
    CR_sem_res_array = np.empty(0) #Array to hold SEMs
    CR_debug_dictionary = {} #Also keep dictionary of averages to help in debugging and cross-condition projection calculations
    analysis_output.write("CROSS-REPLICATE AVERAGES-DISTANCES\n")
    for res_num in sorted(res_num_for_plot):
        average_min_distance = np.average(CR_min_across_reps[res_num])
        sem_min_distance = scipy.stats.sem(CR_min_across_reps[res_num])
        analysis_output.write(f"CR {res_num} average overall: {average_min_distance:.3f}, SEM overall: {sem_min_distance:.3f}\n")
        CR_avg_res_distance_array = np.append(CR_avg_res_distance_array, average_min_distance) #Add new average to the array
        CR_sem_res_array = np.append(CR_sem_res_array, sem_min_distance) #Add new SEM to the array
        CR_debug_dictionary[res_num] = {}
        CR_debug_dictionary[res_num]["avg_distance"] = average_min_distance
        CR_debug_dictionary[res_num]["sem_distance"] = sem_min_distance
    analysis_output.write("\n\n")

    ############
    #Record average projections for each residue
    ############
    CR_avg_res_proj_array = np.empty(0)
    analysis_output.write("CROSS-REPLICATE AVERAGES-CHANNEL PROJ\n")
    for res_num in sorted(res_num_for_plot):
        average_proj = np.average(CR_proj_across_reps[res_num])
        sem_proj = scipy.stats.sem(CR_proj_across_reps[res_num])
        analysis_output.write(f"CR {res_num} average overall: {average_proj:.3f}, SEM overall: {sem_proj:.3f}\n")
        CR_avg_res_proj_array = np.append(CR_avg_res_proj_array, average_proj)
        CR_debug_dictionary[res_num]["avg_proj"] = average_proj
        CR_debug_dictionary[res_num]["sem_proj"] = sem_proj
    analysis_output.write("\n\n")

    #Return arrays of average distance, average projection, and SEM of distance for each residue, as well as number/name translation dictionary and the total frame count and the dictionary of averages/SEMs
    return num_name_map, CR_avg_res_distance_array, CR_avg_res_proj_array, CR_sem_res_array, frame_count_total, CR_debug_dictionary

############
#Run minimum distance analyses and projection analyses for the static structure of interest
############
def static_minimum_distance(simulation_to_analyze, args, analysis_output, dictionary_for_analysis_structures, res_num_for_plot):

    ############
    #Look up the path to the structure to use, record it, and load this structure. Set up the arrays to hold projections and minimum distances
    ############
    structure_to_use = dictionary_for_analysis_structures[simulation_to_analyze]["structure"]
    analysis_output.write(f"STATIC STRUCTURE {simulation_to_analyze}, using {structure_to_use}\n")
    u_structure = MDAnalysis.Universe(structure_to_use, structure_to_use)

    #Set up the minimum distance and projection arrays, to hold the relevant values for each residue selected for analysis
    min_distance_array = np.empty(0)
    proj_array = np.empty(0)

    ############
    #Iterate through the residues of interest and calculate the minimum distance and projection onto the channel axis coordinate for each
    ############
    for r in sorted(res_num_for_plot):

        ############
        #Obtain the minimum distance between all pairs of atoms of this residue on diagonally placed chains
        #Add to the np array for all residues analyzed
        ############
        structure_residue_min_distance_value = minimum_distance_calc(u_structure, r, args, False) 
        min_distance_array = np.append(min_distance_array, structure_residue_min_distance_value)

        ############
        #Obtain the projection of the four relevant CA atoms for this residue onto the channel axis coordinate
        #Add to the np array for all residues analyzed
        ############
        structure_residue_centroid_proj = calc_channel_projection(u_structure, r, args, False)
        proj_array = np.append(proj_array, structure_residue_centroid_proj) #Add on the projection to the growing projection array

        #Record the minimum distance and projection for this residue in the static structure
        analysis_output.write(f"Resid {r} min distance {structure_residue_min_distance_value:.3f}, proj {structure_residue_centroid_proj:.3f}\n")

    #Once minimum distance and projection calculations have been carried out for all residues, return the desired arrays
    return min_distance_array, proj_array

############
#Take in a trajectory directory for a replicate and run MDAnalysis projection and distance calculations, returning the relevant arrays as well as a dictionary mapping residue numbers to residue name/number labels
############
def replicate_distance_and_projection_calc(simulation_selection_name, replicate_directory, args, analysis_output, res_num_for_plot):

    #############
    #Set up- load universe, renumber residues, create data structures
    #############
    #Load the universe for the replicate specified in the function input replicate_directory
    input_xtc = f"{replicate_directory}/md_1.xtc"
    input_tpr = f"{replicate_directory}/md_1.tpr"
    rep_name = replicate_directory.split("/")[-1] #Name for output purposes
    analysis_output.write(f"Values for {rep_name}\n")
    univ_for_analysis = MDAnalysis.Universe(input_tpr, input_xtc)

    #Renumber the residues of the simulation universe so that all four chains have the same residue numbering. The static structure universe does not undergo this processing in static_minimum_distance, because it did not need renumbering in this case, as verified in a Jupyter notebook. Users should check their universes to see whether or not any or all need renumbering. 
    residue_numbers = np.array([at.resid for at in univ_for_analysis.select_atoms("protein")]) #Obtain all residue numbers
    residue_numbers_condensed_list = [r for r, i in groupby(residue_numbers)] #Remove duplicates- this gives 1 resid/residue (https://stackoverflow.com/questions/38065898/how-to-remove-the-adjacent-duplicate-value-in-a-numpy-array)
    residue_numbers_condensed = np.array(residue_numbers_condensed_list)
    residue_numbers_condensed_incremented = residue_numbers_condensed + 1 #Increment protein indices by 1
    residue_numbers_condensed_incremented_mod = residue_numbers_condensed_incremented % 33 #Obtain indices mod 33
    residue_numbers_condensed_incremented_mod[residue_numbers_condensed_incremented_mod == 0] = 33 #Replace 0 with 33 (https://www.w3resource.com/python-exercises/numpy/python-numpy-exercise-88.php, https://stackoverflow.com/questions/3803331/how-can-i-modulo-when-my-numbers-start-from-1-not-zero)
    univ_for_analysis.select_atoms("protein").residues.resids = residue_numbers_condensed_incremented_mod #Set residue IDs to the desired values in the simulation universe

    #Set up dictionary to hold minimum distance and projection data for each residue- keys are "distance" and "proj", and value for each key is a dictionary with keys of residue numbers analyzed and values of the minimum distance (for "distance") and channel axis coordinate projection (for "proj") for each frame considered
    rep_data = {}
    rep_data["distance"] = {}
    rep_data["proj"] = {}

    #Set up dictionary to map residue numbers (keys) to labels of form (resname)(resnumber)- values (this will be used if user wants to label residues on the output plot, in the all_simulation_plot function)
    num_to_name = {}

    #############
    #For each residue the user wants to include in the analysis, obtain the minimum distance array for the simulation and obtain the projection array for the simulation, and also print out the average and SEM for each array
    #############
    for r in sorted(res_num_for_plot):
        
        #Convert the residue number to a label of the form (resname)(resnumber), for the printout and, if wanted, add to the dictionary that will be used for plot labeling
        res_name_for_residue_letter_number_ID = univ_for_analysis.select_atoms(f"protein and resid {r}")[0].resname #This takes the resname for the first atom in the selection corresponding to this residue: because there is one residue number selected, all atoms in the selection should have the same resname
        residue_letter_number_ID = f"{res_name_for_residue_letter_number_ID}{r}"
        if (args.residue_numbers_for_label != None):
            if (r in args.residue_numbers_for_label):
                #If HSD or HSP is in the name, replace with HIS (https://stackoverflow.com/questions/3437059/does-python-have-a-string-contains-substring-method)
                if (("HSD" not in residue_letter_number_ID) and ("HSP" not in residue_letter_number_ID)):
                    residue_letter_number_ID_to_use = residue_letter_number_ID
                if ("HSD" in residue_letter_number_ID):
                    residue_letter_number_ID_to_use = residue_letter_number_ID.replace("HSD", "HIS")
                if ("HSP" in residue_letter_number_ID):
                    residue_letter_number_ID_to_use = residue_letter_number_ID.replace("HSP", "HIS")
                num_to_name[r] = residue_letter_number_ID_to_use
        
        #Run minimum distance function, obtain minimum distance array and times, and store the minimum distance array in rep_data["distance"] with the key corresponding to the residue number
        rep_data["distance"][r], times = minimum_distance_calc(univ_for_analysis, r, args, True)
        analysis_output.write(f"{residue_letter_number_ID_to_use} Average min heavy atom-heavy atom distance is {np.average(rep_data['distance'][r]):.3f} A SEM is {scipy.stats.sem(rep_data['distance'][r]):.3f} A\n") #Record average and SEM

        #Run average placement along channel axis function, obtain the projection array and times, and store the projection array in rep_data["proj"] with the key corresponding to the residue number
        rep_data["proj"][r], times_from_channel = calc_channel_projection(univ_for_analysis, r, args, True)
        analysis_output.write(f"{residue_letter_number_ID_to_use} Average projection is {np.average(rep_data['proj'][r]):.3f} A SEM is {scipy.stats.sem(rep_data['proj'][r]):.3f} A\n") #Record the average and SEM

        #Times should be equal- this is a check to confirm the minimum distance array and the projection array correspond to the same frames
        if (not np.array_equal(times, times_from_channel)): #Use np to confirm the times are all equal (https://docs.scipy.org/doc/numpy/reference/generated/numpy.array_equal.html)
            warnings.warn("***Times indicate different frames pulled for distance and projection calculations")
            analysis_output.write("***Times indicate different frames pulled for distance and projection calculations\n")
            print(times)
            print(times_from_channel)

    #Return the dictionary mapping residue numbers to names, as well as the dictionaries of each residue's minimum distance array and projection array
    return num_to_name, rep_data["distance"], rep_data["proj"]

############
#MDAnalysis call to calculate the minimum distances for a selected residue on diagonally placed chains, and to find and return the minimum of the minima
#Called by minimum_distance_calc for a frame or a static structure
############
def mdanalysis_individual_min_distance_value(selection_1_list, selection_2_list, univ_for_cmd):

    #Calculate distance between all atoms for each selection list
    distance_selection_1 = MDAnalysis.lib.distances.distance_array(univ_for_cmd.select_atoms(selection_1_list[0]).positions, univ_for_cmd.select_atoms(selection_1_list[1]).positions) #As the documentation explains, this command calculates the pairwise distances for all atoms in the two selections (https://www.mdanalysis.org/docs/documentation_pages/lib/distances.html#fast-distance-array-computation-mdanalysis-lib-distances)
    distance_selection_2 = MDAnalysis.lib.distances.distance_array(univ_for_cmd.select_atoms(selection_2_list[0]).positions, univ_for_cmd.select_atoms(selection_2_list[1]).positions)

    #Find the minimum distance between all atoms for each pair of selections, and, of these two minimum values (one for selection_1_list and one for selection_2_list), take the minimum. Then return this minimum-of-minima value
    overall_minimum_distance_value = min(np.min(distance_selection_1), np.min(distance_selection_2))
    return overall_minimum_distance_value

############
#Determine minimum distance (for a structure) or minimum distance time series (for a trajectory) between residues of the same input number (residue_number) on diagonally placed chains
#This can be run for a trajectory (is_trajectory argument set to True, from the function replicate_distance_and_projection_calc) or a static structure (is_trajectory argument set to False, from the function static_minimum_distance)
############
def minimum_distance_calc(univ_to_analyze, residue_number, args, is_trajectory):
    
    ############
    #Select the relevant atoms- heavy atoms for the residue of interest, on chains oriented diagonally to each other. The selection differs slightly for the static structures and the trajectories, given differences in chain naming.
    ############
    if is_trajectory:
        residue_selections_1 = [f"segid seg_0_PROA and resid {residue_number} and not (name H* or name [123]H or type H)", f"segid seg_2_PROC and resid {residue_number} and not (name H* or name [123]H or type H)"] #This syntax enables selecting only the heavy atoms (https://www.mdanalysis.org/docs/documentation_pages/analysis/hbond_analysis.html?highlight=select%20heavy%20atoms)
        residue_selections_2 = [f"segid seg_1_PROB and resid {residue_number} and not (name H* or name [123]H or type H)", f"segid seg_3_PROD and resid {residue_number} and not (name H* or name [123]H or type H)"]

    if not is_trajectory:
        residue_selections_1 = [f"segid PROA and resid {residue_number} and not (name H* or name [123]H or type H)", f"segid PROC and resid {residue_number} and not (name H* or name [123]H or type H)"]
        residue_selections_2 = [f"segid PROB and resid {residue_number} and not (name H* or name [123]H or type H)", f"segid PROD and resid {residue_number} and not (name H* or name [123]H or type H)"]

    ############
    #Using the selection lists, determine minimum distance for the static structure or minimum distance list and associated time list for the trajectory
    ############
    #If the input universe corresponds to a simulation, generate and return minimum distance array and array of associated times
    if is_trajectory:
    
        #Set up lists to hold, for each frame, the minimum distance and the time- to convert to np arrays later
        time_list = [] #Time list
        min_distance_list = [] #Minimum distance list
        
        #Analyzing the desired frames, obtain the minimum pairwise distance for the selections and store the minimum pairwise distance
        for f in univ_to_analyze.trajectory[3000::args.frame_stride]: #Analyze the trajectory starting at 30 ns (consulted https://github.com/MDAnalysis/MDAnalysisCookbook/blob/master/examples/blocks.py, https://www.mdanalysis.org/MDAnalysisTutorial/trajectories.html?highlight=select%20frames has stride)

            #Call the function calculating minimum distances
            min_distance_for_frame = mdanalysis_individual_min_distance_value(residue_selections_1, residue_selections_2, univ_to_analyze)

            #Add the minimum distance calculated to the minimum distance list, and also record the time
            min_distance_list.append(min_distance_for_frame)
            time_list.append(univ_to_analyze.trajectory.time) #Pull time (https://www.mdanalysis.org/MDAnalysisTutorial/trajectories.html?highlight=trajectory%20times)

        #Change lists to arrays after population, for efficiency (https://stackoverflow.com/questions/7332841/add-single-element-to-array-in-numpy)
        min_distance_array = np.array(min_distance_list)
        time_array = np.array(time_list)

        #Return the minimum distance array and time array
        return min_distance_array, time_array

    #If the input universe corresponds to a static structure, generate and return minimum distance value
    if not is_trajectory:
        #Call the function calculating minimum distances
        min_distance_for_structure = mdanalysis_individual_min_distance_value(residue_selections_1, residue_selections_2, univ_to_analyze)
        return min_distance_for_structure

############
#Channel axis projection code written by Michiel Niesen. Many thanks for sharing!
############
def channelaxis(dirv,cv):
    '''
    Determine a normalized vector along the channel axis, based on a direction and a center point.
    '''
    r = np.dot(dirv,cv)
    if r<0:
        return -1*cv/npla.norm(cv)
    else:
        return cv/npla.norm(cv)

############
#Run MDAnalysis call to obtain projections of a specific residue for a specific structure on the channel axis- this will be either for one frame of a trajectory (calc_channel_projection will then include the output in an array for the trajectory) or for a single structure, called by calc_channel_projection function
############
def mdanalysis_individual_projection_value(center, direction, calphas, residue_selection_to_place):
    ############
    #Get the channel axis, many thanks to Michiel Niesen for this code and above channelaxis function!
    ############
    cc = center.centroid() #Center of the channel
    dc = direction.centroid() - cc #Direction of the channel axis vector
    p1, p2, p3 = calphas.principal_axes() #Principal axes of all Calpha atoms
    caxis = channelaxis(dc,p3) #Calculate normalized vector along channel axis

    ############
    #Obtain the centroid of the four relevant CAs and project along the channel axis
    ############
    residue_centroid = residue_selection_to_place.centroid() - cc # Coordinate of this centroid- 4 CA atoms- relative to the center as defined by Ala17 CA atom centroid, thanks so much Michiel!
    residue_centroid_proj = np.dot(residue_centroid, caxis) # Projection on the channel axis, thanks so much Michiel!

    return residue_centroid_proj

############
#Set up MDAnalysis call to obtain projections of a residue on the channel axis
#Many many thanks to Michiel Niesen for sharing!
#This can be run for a trajectory (is_trajectory argument set to True, from replicate_distance_and_projection_calc function) or a static structure (is_trajectory argument set to False, from static_minimum_distance function)
############
def calc_channel_projection(univ_to_analyze_proj, residue_number, args, is_trajectory):

    ############
    #Set up group selections, many thanks to Michiel Niesen for this code
    ############
    center = univ_to_analyze_proj.select_atoms('protein and resname ALA and resid 17 and name CA') #ALA 17 Calphas
    direction = univ_to_analyze_proj.select_atoms('protein and resname TRP and resid 23 and name CA') #TRP 23 Calphas
    calphas = univ_to_analyze_proj.select_atoms('protein and name CA') #All Calphas
    residue_selection_to_place = univ_to_analyze_proj.select_atoms(f"protein and resid {residue_number} and name CA") #Calphas of the residue of interest

    ############
    #If a trajectory is being analyzed, obtain projections for each time desired (starting at 30 ns) by iterating through the trajectory
    ############
    if is_trajectory:

        time_list = [] #List for times
        proj_list = [] #List for projections

        #Iterate through trajectory- at each point obtain the location of the relevant residue with respect to the channel axis
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
#Generate the plot of projected channel axis coordinate versus minimum distance (for trajectories and structures), add lines and labels noting selected residues' channel axis coordinates if wanted, and format and save the plot
############
def all_simulation_plot(args, analysis_output, all_prot_state_dictionary, dictionary_of_colors, res_numbers_for_plot, num_name_linking_dictionary, static_dictionary, static_info, projection_dictionary):

    ############
    #Plot the projected channel axis coordinate against the minimum distance for each residue analyzed for each simulation condition requested (relevant trajectory and structure)
    ############
    for ps in sorted(all_prot_state_dictionary.keys()): #Iterate over simulation conditions analyzed, adding data for each one to the plot

        #Plot minimum distances for the trajectory (https://matplotlib.org/api/_as_gen/matplotlib.pyplot.errorbar.html)
        plt.errorbar(all_prot_state_dictionary[ps]["Minimum_Distances"], all_prot_state_dictionary[ps]["Projections"], xerr = all_prot_state_dictionary[ps]["SEM"], yerr = None, markersize = 10, color = dictionary_of_colors[ps]["data"], label = ps, ecolor = dictionary_of_colors[ps]["uncertainty"])

        #Record trajectory data that was plotted
        analysis_output.write(f"{ps} DATA")
        analysis_output.write("\n")
        for xv, yv, errv in zip(all_prot_state_dictionary[ps]["Minimum_Distances"], all_prot_state_dictionary[ps]["Projections"], all_prot_state_dictionary[ps]["SEM"]):
            analysis_output.write(f"Min distance {xv}, projection {yv}, SEM {errv}")
            analysis_output.write("\n")

        #Plot minimum distances for the structure. No error bars are used here since there is one structure, so there is no uncertainty (https://matplotlib.org/api/_as_gen/matplotlib.pyplot.errorbar.html). The S_ prefix will denote static structures in the legend.
        plt.errorbar(static_dictionary[ps]["Minimum_Distances"], static_dictionary[ps]["Projections"], xerr = None, yerr = None, markersize = 10, color = static_info[ps]["color"], label = f"S_{ps}")

        #Record static structure data that was plotted
        for static_x, static_y in zip(static_dictionary[ps]["Minimum_Distances"], static_dictionary[ps]["Projections"]):
            analysis_output.write(f"Static Min distance {static_x}, projection {static_y}")
            analysis_output.write("\n")

    ############
    #Label axes, set axis label and tick font size, add a legend, and set axis limits if desired
    ############
    #Label axes, set ticks
    plt.ylabel("Channel Axis\nPlacement ($\AA$)", fontsize = 20)
    plt.xlabel("Minimum Distance ($\AA$)", fontsize = 24) #Consulted https://idlangstrom.wordpress.com/2014/12/05/angstrom-and-other-astronomical-symbols-in-python/ for angstrom formatting
    plt.yticks(fontsize = 20) #Change axis label size (https://www.science-emergence.com/Articles/How-to-change-the-size-of-axis-labels-in-matplotlib-/)
    plt.xticks(fontsize = 16) #Change axis label size (https://www.science-emergence.com/Articles/How-to-change-the-size-of-axis-labels-in-matplotlib-/)

    #Add legend, consulted https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot/27355247 for legend help
    plt.legend(loc = "upper left", bbox_to_anchor = (1.02, 1.00)) #Consulted https://matplotlib.org/tutorials/intermediate/legend_guide.html, https://stackoverflow.com/questions/9651092/my-matplotlib-pyplot-legend-is-being-cut-off/42303455
    plt.tight_layout() #Set layout (https://matplotlib.org/3.2.1/tutorials/intermediate/tight_layout_guide.html)

    #Set axes limits if wanted (https://matplotlib.org/api/_as_gen/matplotlib.pyplot.ylim.html?highlight=ylim#matplotlib.pyplot.ylim)
    if (args.minimum_y_for_plot != None):
        plt.ylim(bottom = args.minimum_y_for_plot) 
    if (args.maximum_y_for_plot != None):
        plt.ylim(top = args.maximum_y_for_plot) 
    if (args.minimum_x_for_plot != None):
        plt.xlim(left = args.minimum_x_for_plot)
    if (args.maximum_x_for_plot != None):
        plt.xlim(right = args.maximum_x_for_plot)

    ############
    #If user wants to label certain residues- obtain the projected average coordinate value, draw a line, and label. Note the projection averaged across simulation condition-level averages is drawn (revised from earlier).
    ############
    if (args.residue_numbers_for_label != None):

        minimum_x_for_figure, maximum_x_for_figure = plt.xlim() #Obtain x limits (https://matplotlib.org/3.2.2/api/_as_gen/matplotlib.pyplot.xlim.html?highlight=xlim#matplotlib.pyplot.xlim)
        res_color_for_labeling = (0.60, 0.60, 0.60)

        #Draw a line for each residue desired corresponding to the projection
        for res_num_to_mark in args.residue_numbers_for_label:

            #Obtain the projected-onto-the-z-axis value in the projections dictionary for this residue of interest
            projected_coord_of_res = projection_dictionary[res_num_to_mark]

            #Draw a horizontal line at the channel axis coordinate corresponding to this projection (https://matplotlib.org/3.2.2/api/_as_gen/matplotlib.pyplot.axhline.html?highlight=axhline#matplotlib.pyplot.axhline), and have the line width be small (https:/matplotlib.org/3.1.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D.set_linewidth)
            plt.axhline(projected_coord_of_res, linestyle = "--", linewidth = 1, color = res_color_for_labeling)

            #Write the residue number and associated projected channel axis coordinate value to analysis_output
            analysis_output.write(f"PROJECTION for {res_num_to_mark} {projected_coord_of_res}\n")

            #Add text to associate the appropriate residue name/number identifier with the projection line just drawn (offsets added to maximum_x_for_figure and projected_coord_of_res for cleanliness), using num_name_linking_dictionary to find the correct residue name/number identifier
            plt.text((maximum_x_for_figure - 4), (projected_coord_of_res + 0.50), num_name_linking_dictionary[res_num_to_mark], size = 11, ha = "left", color = res_color_for_labeling) #Add text (https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/arrow_demo.html, https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.text.html)

    ############
    #Assign the plot a title and save it (https://stackoverflow.com/questions/9651092/my-matplotlib-pyplot-legend-is-being-cut-off/42303455)
    ############
    plt.title(f"Minimum Distances Proj. for Residues\nJob {args.title_of_job}", fontsize = 20)
    plt.savefig(f"{args.directory_for_output}/{args.title_of_job}_Minimum_Distance_Proj_Plot.png", bbox_inches = "tight")
    plt.savefig(f"{args.directory_for_output}/{args.title_of_job}_Minimum_Distance_Proj_Plot.eps", bbox_inches = "tight") #Note a message "The PostScript backend does not support transparency; partially transparent artists will be rendered opaque." could appear, workarounds can be found here: https://stackoverflow.com/questions/19638773/matplotlib-plots-lose-transparency-when-saving-as-ps-eps
    plt.savefig(f"{args.directory_for_output}/{args.title_of_job}_Minimum_Distance_Proj_Plot.pdf", bbox_inches = "tight") 


if __name__ == "__main__":

    #############
    #Parse arguments, make output directory/file, set up numbers to plot
    #############
    #Parse arguments
    parser = argparse.ArgumentParser(description = "Plots minimum distance between heavy atoms in opposite chains")
    parser.add_argument("--minimum_residue_number", type = int, help = "Minimum residue number to plot, if a continuous range is wanted. If a continuous range is wanted, this and --maximum_residue_number should both be specified, and -selected_residue_numbers should not be specified.")
    parser.add_argument("--maximum_residue_number", type = int, help = "Maximum residue number to plot, if a continuous range is wanted. If a continuous range is wanted, this and --minimum_residue_number should both be specified, and -selected_residue_numbers should not be specified.")
    parser.add_argument("--selected_residue_numbers", nargs = "+", type = int, help = "Residue numbers to plot, if a continuous range is not wanted. If a continuous range is not wanted, --selected_residue_numbers should be specified and neither --minimum_residue_number nor --maximum_residue_number should be specified.")
    parser.add_argument("--residue_numbers_for_label", nargs = "+", type = int, help = "Residue numbers for labelling on plot, if wanted")
    parser.add_argument("--directory_for_output", required = True, type = str, help = "Directory to be created to hold the results from this analysis")
    parser.add_argument("--title_of_job", required = True, type = str, help = "Name of job. Best to make distinct for each run.")
    parser.add_argument("--dictionary_of_simulations", required = True, type = str, help = "Path of a file (json format) containing a dictionary for accessing simulation data. Keys are shorthand for simulations. Values are dictionaries themselves- named key of 'parent_directory' corresponds to value of directory stem with .tpr and .xtc files.")
    parser.add_argument("--selected_simulations", required = True, nargs = "+", type = str, help = "Keys in dictionary_of_simulations json to be used")
    parser.add_argument("--dictionary_of_colors", required = True, type = str, help = "Path of a file (json format) containing a dictionary of colors for plotting. Keys are those in dictionary_of_simulations. Values are dictionaries themselves- named keys of 'data' and 'uncertainty' correspond each to a list of three floating point numbers representing the colors to use in plotting the values for this set of jobs and the error bars for this set of jobs, respectively.")
    parser.add_argument("--frame_stride", required = True, type = int, help = "Frame stride. This should be small enough that at least two frames are included for each trajectory.")
    parser.add_argument("--minimum_y_for_plot", type = float, help = "Minimum y value for the plot")
    parser.add_argument("--maximum_y_for_plot", type = float, help = "Maximum y value for the plot")
    parser.add_argument("--minimum_x_for_plot", type = float, help = "Minimum x value for the plot")
    parser.add_argument("--maximum_x_for_plot", type = float, help = "Maximum x value for the plot")
    parser.add_argument("--dictionary_of_structures", required = True, type = str, help = "Path of a file (json format) containing a dictionary of structures to be used in minimum distance calculations. Keys are those in --dictionary_of_simulations (or, at a minimum, those requested in --selected_simulations). Values are dictionaries themselves- with key-value pairs of: (1) key 'structure' and value of a path to the structure relevant to that set of simulations and (2) key 'color' and value of a list of three floating point numbers representing the color to use when the minimum distances for this structure are plotted. No uncertainty color is needed because there is only one structure, so there is not uncertainty in this case.")
    args = parser.parse_args()

    #Create output directory and results file (directory cannot already exist: https://www.geeksforgeeks.org/create-a-directory-in-python/)
    os.mkdir(args.directory_for_output)
    analysis_output = open(f"{args.directory_for_output}/Analysis_Output_{args.title_of_job}.txt", "w")

    #Load color dictionary
    with open(args.dictionary_of_colors) as col_dictionary_file_name:
        dictionary_of_colors_for_plotting = json.load(col_dictionary_file_name)

    #Load the directories for analysis dictionary
    with open(args.dictionary_of_simulations) as information_dictionary_file_name:
        directories_for_analysis = json.load(information_dictionary_file_name)

    #Load the static structures for analysis dictionary
    with open(args.dictionary_of_structures) as structure_information_dictionary_file_name:
        dictionary_for_analysis_structures = json.load(structure_information_dictionary_file_name)

    #Record the residue numbers for which the minimum distances on diagonally placed chains will be calculated
    #Either analyze a continuous range of residue numbers (then code in if args.minimum_residue_number != None will be executed) or specify only those wanted (then code in else will be executed)
    #The argparse section above has more details on which arguments should and should not be specified if a continuous range of residue numbers is or is not wanted
    if (args.minimum_residue_number != None):
        res_numbers_for_plot = range(args.minimum_residue_number, args.maximum_residue_number + 1) #Set up list of residue numbers from user input

    else:
        res_numbers_for_plot = args.selected_residue_numbers

    #############
    #Run minimum distance calculations
    #############
    #Set up the results dictionaries: all_prot_state_dictionary holds simulation data. all_prot_state_dictionary_static holds structure data. Also generate all_frame_array to hold frame counts. In addition, all_prot_state_res_named_dictionary holds averaged projection and distance data with the simulation condition as the key and then residues as the next-level keys, to facilitate access for debugging and projection calculations.
    all_prot_state_dictionary = {}
    all_prot_state_dictionary_static = {}
    all_frame_array = np.empty(0)
    all_prot_state_res_named_dictionary = {}
    
    #For each set of jobs of interest, calculate minimum heavy atom-heavy atom distances and projections in the trajectories and the structures
    for simulation_to_analyze in args.selected_simulations:

        #Set up the dictionary entry for this set of simulations, run the minimum distance calculations and obtain the desired arrays, as well as a count of all frames and a dictionary mapping residue numbers to residue names with numbers appended
        all_prot_state_dictionary[simulation_to_analyze] = {}
        num_to_name_dictionary, CRep_avg_res_distance_array, CRep_avg_res_proj_array, CRep_sem_res_array, count_of_all_frames, CRep_dictionary = minimum_distance_for_selection(simulation_to_analyze, args, analysis_output, directories_for_analysis, res_numbers_for_plot)

        #Update dictionaries for simulations and update frame count array
        all_prot_state_dictionary[simulation_to_analyze]["Minimum_Distances"] = CRep_avg_res_distance_array
        all_prot_state_dictionary[simulation_to_analyze]["Projections"] = CRep_avg_res_proj_array
        all_prot_state_dictionary[simulation_to_analyze]["SEM"] = CRep_sem_res_array
        all_prot_state_res_named_dictionary[simulation_to_analyze] = CRep_dictionary
        all_frame_array = np.append(all_frame_array, count_of_all_frames)

        #Run analysis on static structures, obtain the desired arrays, and enter them into all_prot_state_dictionary_static dictionary
        static_mininimum_distance_array, static_projection_on_channel_axis_array = static_minimum_distance(simulation_to_analyze, args, analysis_output, dictionary_for_analysis_structures, res_numbers_for_plot)
        all_prot_state_dictionary_static[simulation_to_analyze] = {}
        all_prot_state_dictionary_static[simulation_to_analyze]["Minimum_Distances"] = static_mininimum_distance_array
        all_prot_state_dictionary_static[simulation_to_analyze]["Projections"] = static_projection_on_channel_axis_array

    #############
    #Generate the minimum distance plot for all jobs
    #############

    #If projections are wanted, calculate them averaged across simulation conditions
    if (args.residue_numbers_for_label != None):

        analysis_output.write("Projections\n")
 
        #Crete dictionary with keys of residue numbers to label. Values will be averaged projections across simulation conditions.
        proj_averaged = dict.fromkeys(args.residue_numbers_for_label)

        #Iterate over residue numbers to label
        for res_to_avg_proj in args.residue_numbers_for_label:
            analysis_output.write(f"Resid {res_to_avg_proj}\n")

            #Empty np array to hold the average projection for each simulation condition
            simulation_condition_averages = np.empty(0)

            #Fill the np array with the average for each condition (i.e., iterate over conditions)
            for state_analyzed in all_prot_state_dictionary.keys(): 

                #Add average projection for each condition to the np array
                projection_for_cond_averaged = all_prot_state_res_named_dictionary[state_analyzed][res_to_avg_proj]["avg_proj"]
                simulation_condition_averages = np.append(simulation_condition_averages, projection_for_cond_averaged)
                analysis_output.write(f"{state_analyzed} projection {projection_for_cond_averaged}\n")

            #Add projection averaged across conditions to proj_averaged
            averaged_proj_value = np.average(simulation_condition_averages)
            proj_averaged[res_to_avg_proj] = averaged_proj_value
            analysis_output.write(f"Resid {res_to_avg_proj} has across conditions an average projection of {averaged_proj_value} A\n")

    if (args.residue_numbers_for_label == None):
        proj_averaged = None

    with open(f"{args.directory_for_output}/{args.title_of_job}_Averages.json", "w") as file_for_saving_av: #Need to use "w" (https://stackoverflow.com/questions/12309269/how-do-i-write-json-data-to-a-file/37795053)
        json.dump(all_prot_state_res_named_dictionary, file_for_saving_av)

    #After all simulation conditions have been analyzed, check that all of them have the same count of observations. If they do not, print out a warning that the average projections across conditions will differ from the average across all pooled frames
    number_of_distinct_obs_counts = len(np.unique(all_frame_array)) #Obtain count of unique values (https://kite.com/python/answers/how-to-count-frequency-of-unique-values-in-a-numpy-array-in-python)
    if (number_of_distinct_obs_counts > 1):
        warnings.warn("***Please note observation counts are not the same for all simulation conditions. Be aware that the cross-condition residue projection average, which is reported on the plot if projections are requested, will be different from the average obtained by weighting each frame equally.")
        analysis_output.write("***Please note observation counts are not the same for all simulation conditions. Be aware that the cross-condition residue projection average, which is reported on the plot if projections are requested, will be different from the average obtained by weighting each frame equally.\n")

    #Now generate plot
    all_simulation_plot(args, analysis_output, all_prot_state_dictionary, dictionary_of_colors_for_plotting, res_numbers_for_plot, num_to_name_dictionary, all_prot_state_dictionary_static, dictionary_for_analysis_structures, proj_averaged)

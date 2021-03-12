import argparse
import glob
from itertools import groupby
import json
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis
import MDAnalysis.analysis.rms
import numpy as np
import os
import pickle

#############
#For each set of replicates, calculate and plot RMSDs
#Print out cross-replicate averages. If the first frame is the same for all replicates, it will be counted as many times in the average as there are replicates. For low frame counts this may be an issue, but the effect should be negligible for large frame counts, as used in analysis when this script was run for the BM2 project.
#############
def rmsd_calculations(directory_dictionary_for_analysis, analysis_output, dictionary_of_rmsd_selections, color_dictionary_for_plotting, args):

    #Create dictionary to hold the RMSD time series
    #The format is key = name of set of jobs, value = dictionary for each job where: 
    #key = name of replicate, value = dictionary for each selection, where:
    #key = name of selection, value = RMSD output transposed
    overall_rmsd_dictionary = {}

    #############
    #For each set of simulations, identify the relevant directories, calculate RMSD time series for each simulation, plot the results, and report the average across replicates
    #############
    for simulation_selection_name in directory_dictionary_for_analysis.keys():
        overall_rmsd_dictionary[simulation_selection_name] = {}

        #Obtain the individual replicate directories for each directory (https://stackoverflow.com/questions/14506491/php-glob-pattern-match-for-arbitray-number-of-digits) and document in output file
        print(f"On {simulation_selection_name}") #Use f-strings to format (https://realpython.com/python-f-strings/)
        analysis_output.write(f"=================\n{simulation_selection_name}\n")
        print(f"{directory_dictionary_for_analysis[simulation_selection_name]['parent_directory']}")
        replicate_directories_1_digit = glob.glob(f"{directory_dictionary_for_analysis[simulation_selection_name]['parent_directory']}[0-9]") #This is written for a directory system where there is a numerical suffix for each replicate, and there are up to 100 replicates
        replicate_directories_2_digit = glob.glob(f"{directory_dictionary_for_analysis[simulation_selection_name]['parent_directory']}[0-9][0-9]")
        replicate_directories = sorted(replicate_directories_1_digit) + sorted(replicate_directories_2_digit)
        analysis_output.write(f"*************\nAnalysis for {simulation_selection_name}\nParent directory {directory_dictionary_for_analysis[simulation_selection_name]['parent_directory']}\n")

        #Also record the path to the structure that will be the reference in the RMSD calculations
        protein_structure_reference_for_rmsd = directory_dictionary_for_analysis[simulation_selection_name]["protein_structure_reference"]

        #Set up replicate-level averages dictionary- keys are selection shorthand labels, and values will be np arrays of averages
        replicate_level_average_dictionary = dict.fromkeys(dictionary_of_rmsd_selections.keys(), np.empty(0)) #https://www.programiz.com/python-programming/methods/dictionary/fromkeys

        #For each replicate run RMSD calculations and add the individual replicate's output to the overall_rmsd_dictionary
        for rep_d in replicate_directories:
            print(f"On {rep_d}")
            analysis_output.write(f"{rep_d}\n")

            #Generate a key name for this replicate and run RMSD calculations. Put output (a dictionary with keys of selection names and values of MDAnalysis output) inside overall_rmsd_dictionary as a value, with the key being the replicate.
            k_name = f"{simulation_selection_name}_Rep_{rep_d.split('_')[-1]}"
            overall_rmsd_dictionary[simulation_selection_name][k_name] = rmsd_calc(rep_d, analysis_output, dictionary_of_rmsd_selections, args, protein_structure_reference_for_rmsd, k_name)

            #Record average for each replicate for each selection in the appropriate np array in replicate_level_average_dictionary
            for rmsd_selection_to_average in replicate_level_average_dictionary.keys():
                replicate_level_average_dictionary[rmsd_selection_to_average] = np.append(replicate_level_average_dictionary[rmsd_selection_to_average], np.average(overall_rmsd_dictionary[simulation_selection_name][k_name][rmsd_selection_to_average][2])) #Details on MDAnalysis output can be found here: https://www.mdanalysis.org/docs/documentation_pages/analysis/rms.html#MDAnalysis.analysis.rms.rmsd

        #For each MDAnalysis selection of the protein analyzed, plot the time series for each replicate and record the average across replicates
        #If the first frame is the same for all replicates, it will be counted as many times in the average as there are replicates. For low frame counts this may be an issue, but the effect should be negligible for large frame counts, as used in analysis when this script was run for the BM2 project.
        for rmsd_selection_for_plot in dictionary_of_rmsd_selections.keys():
            rmsd_plot(overall_rmsd_dictionary[simulation_selection_name], simulation_selection_name, rmsd_selection_for_plot, color_dictionary_for_plotting, args) #Only pass the relevant portion of overall_rmsd_dictionary into the plotting function, not the entire dictionary (e.g., only pass in P01 unrestrained simulations in a dictionary containing P01 unrestrained, P01 restrained, P4 unrestrained, and P4 restrained simulations)
            analysis_output.write(f"Average RMSD for {simulation_selection_name} {rmsd_selection_for_plot} is {np.average(replicate_level_average_dictionary[rmsd_selection_for_plot]):.4f}\n") #Get cross-replicate average for a desired set of simulations and for a desired selection (formatting information here: https://stackoverflow.com/questions/1995615/how-can-i-format-a-decimal-to-always-show-2-decimal-places)

    return overall_rmsd_dictionary

#############
#Carry out RMSD calculation for all requested selections for a given replicate
#Returns a dictionary with keys of the selection name shorthand (args.dictionary_of_selections keys) and values of the MDAnalysis.analysis.rms.RMSD output for each selection
#############
def rmsd_calc(rep_direc, analysis_output, dictionary_of_selections_of_rmsds, args, protein_structure_reference_for_rmsd, rep_shorthand):

    #Load trajectory for this simulation and structure reference for this set of simulations
    input_xtc = f"{rep_direc}/md_1_nojump.xtc"
    input_tpr = f"{rep_direc}/md_1.tpr"
    u_rmsd = MDAnalysis.Universe(input_tpr, input_xtc)
    analysis_output.write(f"RMSD xtc {input_xtc}\n")
    analysis_output.write(f"RMSD tpr {input_tpr}\n")
    analysis_output.write(f"Structure reference {protein_structure_reference_for_rmsd}\n")
    u_reference = MDAnalysis.Universe(protein_structure_reference_for_rmsd, protein_structure_reference_for_rmsd)
    
    #Renumber the residues of the simulation universe so that all four chains have the same residue numbering. The reference universe did not need renumbering in this case, as verified in a Jupyter notebook. Users should check their universes to see whether or not any or all need renumbering.
    residue_numbers = np.array([at.resid for at in u_rmsd.select_atoms("protein")]) #Obtain all residue numbers
    residue_numbers_condensed_list = [r for r, i in groupby(residue_numbers)] #Remove duplicates- this gives 1 resid/residue (https://stackoverflow.com/questions/38065898/how-to-remove-the-adjacent-duplicate-value-in-a-numpy-array)
    residue_numbers_condensed = np.array(residue_numbers_condensed_list)
    residue_numbers_condensed_incremented = residue_numbers_condensed + 1 #Increment protein indices by 1
    residue_numbers_condensed_incremented_mod = residue_numbers_condensed_incremented % 33 #Obtain indices mod 33
    residue_numbers_condensed_incremented_mod[residue_numbers_condensed_incremented_mod == 0] = 33 #Replace 0 with 33 (https://www.w3resource.com/python-exercises/numpy/python-numpy-exercise-88.php, https://stackoverflow.com/questions/3803331/how-can-i-modulo-when-my-numbers-start-from-1-not-zero)
    u_rmsd.select_atoms("protein").residues.resids = residue_numbers_condensed_incremented_mod #Set residue IDs to the desired values in the simulation universe

    #Set up results dictonary that will be returned. Keys are the selection name shorthand (args.dictionary_of_selections keys) and values are the transposed MDAnalysis.analysis.rms.RMSD output for each selection.
    rmsd_dictionary = {}

    #For each selection, run the RMSD calculation
    for rmsd_selection_label, rmsd_selection_text in dictionary_of_selections_of_rmsds.items():

        #Run RMSD calculation, put in dictionary
        #RMSD calculation information source: https://www.mdanalysis.org/docs/documentation_pages/analysis/rms.html#MDAnalysis.analysis.rms.rmsd
        rmsd_to_calc = MDAnalysis.analysis.rms.RMSD(u_rmsd, reference = u_reference, select = rmsd_selection_text, filename = f"{args.directory_for_output}/{rep_shorthand}_{rmsd_selection_label}_RMSD.dat")
        rmsd_to_calc.run(step = args.frame_stride)
        rmsd_dictionary[rmsd_selection_label] = rmsd_to_calc.rmsd.T #Documentation suggested taking the transpose to facilitate figure generation

    #Return the dictionary, to be added to overall_rmsd_dictionary in rmsd_calculations
    return rmsd_dictionary

#############
#Generate plots for the input RMSD selection, for one set of replicates. This is run once for each replicate for each RMSD selection by a loop in the rmsd_calculations function, so RMSD selections are not iterated over in this function.
#############
def rmsd_plot(dictionary_for_plot, analysis_selection_name, rmsd_selection_for_plot, color_dictionary_for_plotting, args):
    #Obtain relevant color list: full job name from args.dictionary_of_simulations is checked against args.dictionary_of_colors keys which are abbreviated
    color_list_key_to_use = [col_key for col_key in color_dictionary_for_plotting.keys() if col_key in analysis_selection_name][0] #Obtain the color dictionary key to use from checking color dictionary keys against the args.dictionary_of_simulations selection names. (The argparse section in the main portion of this script has more details on required key formatting, which enables this code line to work. Though a list is returned, it should only have one item as long as key formatting is done correctly, which is why [0] is used. Substring search information can be found here: https://stackoverflow.com/questions/3437059/does-python-have-a-string-contains-substring-method.)
    color_list_to_use = color_dictionary_for_plotting[color_list_key_to_use] #Obtain the list of colors to use from the key found to be the appropriate one

    #Iterate through the keys, each corresponding to a replicate, in the dictionary. For each, plot the times on the x axis and RMSD values on the y axis. The use of enumerate enables having a different color for each replicate, since it accesses the different colors listed in color_list_to_use. Please see the notes in the argparse section in the main portion of this script regarding the length of color_list_to_use.
    #Time and RMSD are extracted from the MDAnalysis.analysis.rms.RMSD output following the documentation/examples (https://www.mdanalysis.org/docs/documentation_pages/analysis/rms.html#MDAnalysis.analysis.rms.rmsd).
    for index_for_color, simulation_label in enumerate(sorted(dictionary_for_plot.keys())):
        plt.plot(dictionary_for_plot[simulation_label][rmsd_selection_for_plot][1] / 1000, dictionary_for_plot[simulation_label][rmsd_selection_for_plot][2], color = color_list_to_use[index_for_color], label = simulation_label, alpha = 0.9)

    #Set axis limits if desired, format axis/label sizes, title the plot, add a legend, and save
    #Legend resources for this/other scripts: https://stackoverflow.com/questions/9651092/my-matplotlib-pyplot-legend-is-being-cut-off/42303455, https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot/27355247, https://matplotlib.org/tutorials/intermediate/legend_guide.html
    #Set axis limits if desired (https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.ylim.html)
    if (args.minimum_y_for_plot != None):
        plt.ylim(bottom = args.minimum_y_for_plot)
    if (args.maximum_y_for_plot != None):
        plt.ylim(top = args.maximum_y_for_plot)
    if (args.minimum_x_for_plot != None):
        plt.xlim(left = args.minimum_x_for_plot)
    if (args.maximum_x_for_plot != None):
        plt.xlim(right = args.maximum_x_for_plot)
    plt.title(f"RMSD Time Series St {args.frame_stride}\n{analysis_selection_name} Sim.s\n{rmsd_selection_for_plot} Atoms\n{ args.title_of_job}", fontsize = 18)
    plt.xlabel("Time(ns)", fontsize = 22) #Times initially were in ps (https://www.mdanalysis.org/MDAnalysisTutorial/trajectories.html?highlight=trajectory%20times)
    plt.ylabel("RMSD ($\AA$)", fontsize = 22) #Consulted https://idlangstrom.wordpress.com/2014/12/05/angstrom-and-other-astronomical-symbols-in-python/ for angstrom formatting
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 22)
    plt.legend(loc = "upper left", bbox_to_anchor = (1.02, 1.00), prop={"size": 8}) #Place legend (https://stackoverflow.com/questions/7125009/how-to-change-legend-size-with-matplotlib-pyplot)
    plt.tight_layout() #Set layout (https://matplotlib.org/3.2.1/tutorials/intermediate/tight_layout_guide.html)
    plt.savefig(f"{args.directory_for_output}/{analysis_selection_name}_{rmsd_selection_for_plot}_RMSD_{args.title_of_job}.png")
    plt.savefig(f"{args.directory_for_output}/{analysis_selection_name}_{rmsd_selection_for_plot}_RMSD_{args.title_of_job}.pdf")
    plt.savefig(f"{args.directory_for_output}/{analysis_selection_name}_{rmsd_selection_for_plot}_RMSD_{args.title_of_job}.eps") #Save figure, note a message "The PostScript backend does not support transparency; partially transparent artists will be rendered opaque." could appear, workarounds can be found here: https://stackoverflow.com/questions/19638773/matplotlib-plots-lose-transparency-when-saving-as-ps-eps

    #Clear/close in preparation for next plot (https://stackoverflow.com/questions/17106288/matplotlib-pyplot-will-not-forget-previous-plots-how-can-i-flush-refresh)
    plt.clf()
    plt.cla()
    plt.close()

if __name__ == "__main__":

    #############
    #Parse arguments, make output directory/file, load relevant dictionaries
    #############
    #Parse arguments
    parser = argparse.ArgumentParser(description = "Plots RMSD time series")
    parser.add_argument("--directory_for_output", required = True, type = str, help = "Directory to be created to hold the results from this analysis")
    parser.add_argument("--title_of_job", required = True, type = str, help = "Name of job. Best to make distinct for each run.")
    parser.add_argument("--dictionary_of_simulations", required = True, type = str, help = "Path of a file (json format) containing a dictionary for accessing simulation data. Keys are shorthand for simulations. Values are dictionaries themselves- named keys of 'parent_directory' and 'protein_structure_reference' correspond to values of directory stem with .tpr and .xtc files and the structure reference for the RMSD calculations, respectively.")
    parser.add_argument("--dictionary_of_colors", required = True, type = str, help = "Path of a file (json format) containing a dictionary of colors for plotting. Keys are shorthand for simulations (abbreviated or full versions of the keys in dictionary_of_simulations), and values are lists of colors (each color is expressed as a list of three floating point numbers, and there should be at least as many colors given as there are replicates for the relevant set(s) of simulations for which the list will be accessed in the rmsd_plot function). Please note that each --dictionary_of_simulations key must have or have as a substring exactly one of the --dictionary_of_colors keys, in order for color selections to be accessed correctly in the rmsd_plot function. For instance, --dictionary_of_simulations keys could be 'P01_Restrained', 'P01_No_Restraints', 'P4_Restrained', and 'P4_No_Restraints', and dictionary_of_colors keys could be 'P01' and 'P4'")
    parser.add_argument("--dictionary_of_selections", required = True, type = str, help = "Path of a file (json format) containing a dictionary of selections. Keys are shorthands for selections, and values are the MDAnalysis-compatible selections.")
    parser.add_argument("--frame_stride", required = True, type = int, help = "Value for striding")
    parser.add_argument("--minimum_y_for_plot", type = float, help = "Minimum y value for plot")
    parser.add_argument("--maximum_y_for_plot", type = float, help = "Maximum y value for plot")
    parser.add_argument("--minimum_x_for_plot", type = float, help = "Minimum x value for plot")
    parser.add_argument("--maximum_x_for_plot", type = float, help = "Maximum x value for plot")
    args = parser.parse_args()

    #Create output directory and results file (directory cannot already exist: https://www.geeksforgeeks.org/create-a-directory-in-python/)
    os.mkdir(args.directory_for_output)
    analysis_output = open(f"{args.directory_for_output}/Analysis_Output_{args.title_of_job}.txt", "w")

    #Load color, selections, directory dictionaries
    with open(args.dictionary_of_colors) as col_dictionary_file_name:
        color_dictionary_for_plot = json.load(col_dictionary_file_name)

    with open(args.dictionary_of_selections) as selections_of_atoms_dictionary_file_name:
        selections_dictionary = json.load(selections_of_atoms_dictionary_file_name)

    with open(args.dictionary_of_simulations) as simulations_dictionary_file_name:
        simulations_dictionary = json.load(simulations_dictionary_file_name)

    #############
    #Run the RMSD calculations and plot the time series.
    #############
    time_series_dictionary = rmsd_calculations(simulations_dictionary, analysis_output, selections_dictionary, color_dictionary_for_plot, args)

    #Save out plot data with "wb" setting due to error (https://stackoverflow.com/questions/5512811/builtins-typeerror-must-be-str-not-bytes)
    with open(f"{args.directory_for_output}/{args.title_of_job}_RMSD_Analysis_Overall_Dictionary.json", "wb") as file_for_saving: pickle.dump(time_series_dictionary, file_for_saving)

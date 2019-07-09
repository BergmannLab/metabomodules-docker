# Copyright 2019 Bergmann's Lab UNIL <mattia.tomasoni@unil.ch>
#
#    This file is part of METABOMODULES Tool.
#
#    METABOMODULES Tool is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    METABOMODULES Tool is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    GNU General Public License: <https://www.gnu.org/licenses/>.

import datetime
import os
import argparse
import src.pca.pca_pseudospectra as pca_pseudospectra
import src.acp.acp_wrp as acp_wrp
import src.isa.isawrp as isa_wrp
import src.isa.attractor as attractor
import src.filter.filtersig_wrp as filtersig_wrp
import metabomodules_config

octave_root_dir = metabomodules_config.resources['octave_root_dir'] + '/'
ISApaper_root_dir = metabomodules_config.resources['current_dir'] + '/'

isa_root_dir = ISApaper_root_dir + 'src/isa/'
stocsy_root_dir = ISApaper_root_dir + 'src/stocsy/ForWrapperCode/'
filter_root_dir = ISApaper_root_dir + 'src/filter/'
metabomatching_root_dir = ISApaper_root_dir + 'src/metabomatching/'



def prepare_output_dir(output_root_dir, redo_flag):
    if output_root_dir[0]!="/": 
        output_root_dir = ISApaper_root_dir + output_root_dir
    if output_root_dir[-1]!="/": 
        output_root_dir = output_root_dir + "/"
    
    if redo_flag==False:
        output_dir = output_root_dir + datetime.datetime.now().strftime("%Y-%m-%d__%H%M%S") + "/"
        os.system("mkdir -p " + output_dir)
    else:
        previous_outputs = os.listdir(output_root_dir);
        output_dir = output_root_dir + previous_outputs[-1] + "/";
    return output_dir
    
def generate_parameters_file(parameters_file, number_metabomatching_permutations,
                             numb_samples = None, crScale = None):
    with open(parameters_file, "w") as file:
        file.write("reference\tumdb\n")
        file.write("dsh\t0.025\n")
        file.write("significant\t4\n")
        file.write("suggestive\t3\n")
        file.write("nshow\t16\n")
        file.write("variant\tpm\n")
        file.write("n_permutation\t" + str(number_metabomatching_permutations) + "\n")
        if numb_samples != None:
            file.write("samplesize\t" + str(numb_samples) + "\n")
        if crScale != None:
            file.write("crscale\t" + crScale + "\n")

def run_PCA(dataset_name, input_file, output_dir, number_pseudospectra, number_metabomatching_permutations):
    pca_output_dir = output_dir + 'ps.pca.' + dataset_name + '/'
    os.mkdir(pca_output_dir)
    pca_pseudospectum = pca_output_dir + 'pca.' + dataset_name + '.pseudospectrum.tsv'
    pca_pseudospectra.generate_pca_pseudospectra(input_file, pca_pseudospectum, number_pseudospectra)
    pca_parameters_file = pca_output_dir + 'parameters.in.tsv'
    generate_parameters_file(pca_parameters_file, number_metabomatching_permutations)
    
def run_ACP(dataset_name, input_file, output_dir, number_pseudospectra, number_metabomatching_permutations, OffDiagDist, remNeigbPairsFlag):
    stocsy_output_dir = output_dir + 'ps.acp.' + dataset_name + '/'
    os.mkdir(stocsy_output_dir)
    (numb_samples, crScale) = acp_wrp.main(input_file, stocsy_output_dir, number_pseudospectra, OffDiagDist, remNeigbPairsFlag)
    stocsy_parameters_file = stocsy_output_dir + 'parameters.in.tsv'
    generate_parameters_file(stocsy_parameters_file, number_metabomatching_permutations, numb_samples, crScale)

def run_ISA(dataset_name, input_file, output_dir, number_metabomatching_permutations, args):
    isa_output_dir = output_dir + 'ps.isa.' + dataset_name + '/'
    os.mkdir(isa_output_dir)
    clean_previous_run = "rm -rf " + isa_root_dir + "*.tsv "+ isa_root_dir + "*.csv"
    os.system(clean_previous_run)
    move_to_isa_dir = "cp " + input_file + " " + isa_root_dir
    os.system(move_to_isa_dir)
    os.chdir(isa_root_dir)
    
    # initial run
    args.inpfile = dataset_name + ".csv"; args.outfile = dataset_name + ".csv"
    args.header = True; args.label = True; args.gopseudo = True; args.seedsparsity = 3; args.nt = True
    args.thc = "1,2,3,4,5,6,7"; args.thr = "1,2,3,4,5,6,7"
    args.sgc = 0; args.sgr = 1; args.dconv = 0.99; args.dsame = 0.50; args.nseed = 250
    isa_wrp.main(input_file, args)
    # turn off sweeping and run iterations
    for i in range(11, 31):
        args.outfile = dataset_name + ".attract" + str(i) + ".csv"
        args.dsame = 0.8; args.nseed = 500; args.nosweep = True; args.nopurge = True;
        isa_wrp.main(input_file, args)
    attractor.main(dataset_name + ".colscore.tsv")

    move_back = "cp " + isa_root_dir + "isa.*pseudospectrum.tsv " + isa_root_dir + "isa.*info.tsv " + isa_output_dir
    os.system(move_back)
    # create metabomatching parameters file
    isa_parameters_file = isa_output_dir + 'parameters.in.tsv'
    generate_parameters_file(isa_parameters_file, number_metabomatching_permutations)

def run_metabomatching(output_dir, use_octave, number_p):
    clean_previous_run = "rm -rf " + metabomatching_root_dir + "ps.* "
    move_to_metabomatching_dir = "mv " + output_dir + "ps.* " + metabomatching_root_dir
    enter_metabomatching_dir = "cd " + metabomatching_root_dir
    move_back = "mv " + metabomatching_root_dir + "ps.* " + output_dir
    if use_octave:
        print("Run metabomatching with Octave")
        run_metabomatching = octave_root_dir + "octave metabomatching.m"
    else:
        print("Run metabomatching with MATLAB") 
        run_metabomatching = matlab_root_dir + "matlab -nodisplay -nosplash -nodesktop -r \"run('" \
            + metabomatching_root_dir + "metabomatching.m');quit\""
    print("NUMBER OF PERMUTATIONS: " + str(number_p))
    os.system(clean_previous_run)
    os.system(move_to_metabomatching_dir)
    os.system(enter_metabomatching_dir + " && " + run_metabomatching)
    os.system(move_back)

def run_filter(output_dir, z_score_threshold, adj_score_threshold, redo_flag):
    if redo_flag: # only re-run the filter on the output of the previous run
        output_dir_previous_run = output_dir
        output_dir = output_dir + "original/" # point to the un-filtered output of the previous run
        move_back_filtered = "mv " + filter_root_dir + "z_score* " + output_dir_previous_run # filtered results
        move_back_original = "mv " + filter_root_dir + "ps* " + output_dir # un-filtered, original, results
        move_back = move_back_filtered + " && " + move_back_original
    else: # normal filter run
        move_back = "mv " + filter_root_dir + "z_score* " + filter_root_dir + "original " + output_dir
    
    clean_previous_run = "rm -rf " + filter_root_dir + "ps.* " + filter_root_dir + "z_score* " + filter_root_dir + "original"
    os.system(clean_previous_run)
    move_to_filter_dir = "mv " + output_dir + "ps.* " + filter_root_dir
    os.system(move_to_filter_dir)
    filtered_folder = filtersig_wrp.main(filter_root_dir, z_score_threshold, adj_score_threshold, redo_flag)
    if redo_flag: # in case the user re-run the filter with the same thresholds
        remove_duplicates = "rm -rf " + output_dir_previous_run + filtered_folder
        os.system(remove_duplicates) 
    os.system(move_back)
    return filtered_folder + "/"



def generate_pseudospectra(output_dir, number_pseudospectra, number_metabomatching_permutations, args):
    for dataset in args.datasets:
        dataset_name = dataset.split("/")[-1][0:-4]
        input_file = "___"
        if dataset[0]!="/": 
            input_file = ISApaper_root_dir + dataset
        else:
            input_file = dataset
        #print("\n---------------------------RUNNING ACP-------------------------------\n")
        run_ACP(dataset_name, input_file, output_dir, 
                   number_pseudospectra, number_metabomatching_permutations, args.OffDiagDist, args.remNeigbPairsFlag)
        #print("\n---------------------------RUNNING PCA-------------------------------\n")
        #run_PCA(dataset_name, input_file, output_dir, 
        #        number_pseudospectra, number_metabomatching_permutations)
        #print("\n---------------------------RUNNING ISA-------------------------------\n")
        #run_ISA(dataset_name, input_file, output_dir, number_metabomatching_permutations, args)

def run_pipeline(args):
    number_metabomatching_permutations = 99
    output_dir = prepare_output_dir(args.output_root_dir, args.redo_flag)
    if args.redo_flag==False:
        # RUN ALL METHODS ON ALL DATASETS
        print("\n------RUNNING MODULARIZATION METHODS: GENERATING PSEUDOSPACTRA-------\n")
        generate_pseudospectra(output_dir, args.number_pseudospectra, number_metabomatching_permutations, args)
        # RUN METABOMATCHING ON ALL PSEUDOSPECTRA
        #print("\n-------------------METABOMATCHING INITIAL RESULTS--------------------\n")
        #run_metabomatching(output_dir, args.use_octave, number_metabomatching_permutations)
        # RUN FILTER
        #print("\n------------------------FILTERING RESULTS----------------------------\n")
        #filtered_folder = run_filter(output_dir, args.z_score_threshold, args.adj_score_threshold, args.redo_flag)
        #output_dir = output_dir + filtered_folder
        # RE-RUN METABOMATCHING ON FILTERED PSEUDOSPECTRA
        #print("\n------------------METABOMATCHING FILTERED RESULTS--------------------\n")
        #run_metabomatching(output_dir, args.use_octave, number_metabomatching_permutations)
    else:
        # RUN-FILTER
        print("\n-----------------------UN-FILTERING RESULTS---------------------------\n")
        filtered_folder = run_filter(output_dir, args.z_score_threshold, args.adj_score_threshold, args.redo_flag)
        output_dir = output_dir + filtered_folder
        # RE-RUN METABOMATCHING ON FILTERED PSEUDOSPECTRA
        print("\n------------------METABOMATCHING FILTERED RESULTS--------------------\n")
        run_metabomatching(output_dir, args.use_octave, number_metabomatching_permutations)

if __name__ == "__main__":
    main_parser = argparse.ArgumentParser(description='METABOMODULES')

    # GENERAL PARAMETERS
    main_parser.add_argument('-N', action="store", dest="number_pseudospectra", type=int, 
                    help='Number of pseudospectra')
    group = main_parser.add_mutually_exclusive_group()
    group.add_argument('--octave', action='store_true', default=True, dest='use_octave',
                    help='Use octave to run metabomatching')
    group.add_argument('--matlab', action='store_false', default=True, dest='use_octave',
                    help='Use matlab to run metabomatching')
    main_parser.add_argument('--output', action="store", dest="output_root_dir",
                    help='Output folder')
    required_arguments = main_parser.add_argument_group('required named arguments')
    required_arguments.add_argument('--input', action='append', dest='datasets', default=[],
                    help='Add file to the list of datasers be processed', required=True)
    # ACP PARAMETERS
    main_parser.add_argument('--ACP_remNeigbPairsFlag', action="store_true", default=True, dest="remNeigbPairsFlag",
                    help='ACP parameter: remove pairs that lie in the same neighborhood as other pairs')
    main_parser.add_argument('--ACP_OffDiagDist', action="store", default=0.1, dest="OffDiagDist", type=float, 
                    help='Off diagonal distance for correlation matrix')
    
    # ISA PARAMETERS
    main_parser.add_argument('--ISA_seedmatrix',dest='seedfile',default=None)
    main_parser.add_argument('--ISA_dsame',dest='dsame',default=0.80)
    main_parser.add_argument('--ISA_dconv',dest='dconv',default=0.975)
    main_parser.add_argument('--ISA_nseed',dest='nseed',default=100)
    main_parser.add_argument('--ISA_seedsparsity',dest='seedsparsity',type=int,default=0)
    main_parser.add_argument('--ISA_maxiter',dest='maxiter',default=50)
    main_parser.add_argument('--ISA_sgc',dest='sgc',default=0)
    main_parser.add_argument('--ISA_sgr',dest='sgr',default=1)
    main_parser.add_argument('--ISA_thc',dest='thc',default='1,2,3')
    main_parser.add_argument('--ISA_thr',dest='thr',default='1,2,3')
    main_parser.add_argument('--ISA_norm',dest='norm',default='double')
    main_parser.add_argument('--ISA_nt',action='store_true',dest='nt',default=False)
    main_parser.add_argument('--ISA_inputhasheaders',action='store_true',dest='header',default=False)
    main_parser.add_argument('--ISA_inputhaslabels',action='store_true',dest='label',default=False)
    main_parser.add_argument('--ISA_nopurge',action='store_true',dest='nopurge',default=False)
    main_parser.add_argument('--ISA_quiet',action='store_true',dest='quiet',default=False)
    main_parser.add_argument('--ISA_nosweep',action='store_true',dest='nosweep',default=False)
    main_parser.add_argument('--ISA_onefile',action='store_true',dest='onefile',default=False)
    main_parser.add_argument('--ISA_gopseudo',action='store_true',dest='gopseudo',default=False)
    
    # FILTER PARAMETERS
    main_parser.add_argument('--FILETR_z_score_threshold', action="store", default=4, dest="z_score_threshold", type=float, 
                    help='pseudospectra with values below this threshold are not processed')
    main_parser.add_argument('--FILETR_adj_score_threshold', action="store", default=2, dest="adj_score_threshold", type=float, 
                    help='pseudospectra that score below this threshold are not processed')
    main_parser.add_argument('--FILETR_redo', action="store_true", default=False, dest="redo_flag",
                    help='redo filtering on the output of a previously processed dataset')
    
    print("\n------------------STARTING METABOMODULES PIPELINE--------------------\n")
    run_pipeline(main_parser.parse_args())
    print("\n------------------FINISHED METABOMODULES PIPELINE--------------------\n")

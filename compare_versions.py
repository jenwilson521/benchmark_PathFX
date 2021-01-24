import numpy as np
import pandas as pd
import sys
import os
import argparse
import multiprocessing as mp

print("\n")

# Setting up argparse
parser = argparse.ArgumentParser(description="Calculate differences between between identified phenotypes and genes "
                                             "to help compare PathFX versions")

command_group_optional_args = parser.add_mutually_exclusive_group()
command_group_optional_args.add_argument('-d', '--drugs', metavar='drug_list', type=str, nargs='+',
                                         help='List of drugs to compare with PathFX versions. Compares '
                                              'genes by default. Use --phenotypes to compare with phenotypes.')
command_group_optional_args.add_argument('-r', '--random-drugs', metavar="n_drugs", type=int,
                                         help='Number of drugs to randomly compare')
parser.add_argument('-p', '--phenotypes', help='Compares PathFX using the amount of phenotypes it identifies',
                    action='store_true')

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


# Returns how many more genes v2 identified than v1
def compare_one_drug(v1db, v2db, print_stuff=True):
    # Ask until valid input
    while True:

        # Check if text file
        if v1db[-4:] == ".txt":
            break
        else:
            print(
                "\nPlease enter the path to a .txt file or drag and drop it into the terminal window for the first "
                "drug.")

        # Request the v1 merged_neighborhood__assoc_database_sources data for the drug
        v1db = input("Path to v1 <drug>_merged_neighborhood__assoc_database_sources_.txt: \n")

        # Strip extra spaces if they exist
        v1db = v1db.rstrip()

    while True:

        # Check if text file
        if v2db[-4:] == ".txt":
            break
        else:
            print(
                "\nPlease enter the path to a .txt file or drag and drop it into the terminal window for the second "
                "drug.")

        # Request the v2 merged_neighborhood__assoc_database_sources data for the drug
        v2db = input("\nPath to v2 <drug>_merged_neighborhood__assoc_database_sources_.txt: \n")

        # Strip extra spaces if they exist
        v2db = v2db.rstrip()

    data_frame_v1 = pd.read_csv(v1db, delimiter="\t", index_col=False)

    # v2 algorithm does not give the phenotype column a header, so headers have to be overridden
    data_frame_v2 = pd.read_csv(v2db, delimiter="\t", index_col=False, header=0, names=['Gene', 'CUI', 'Phenotype',
                                                                                        'Source Databases'])

    # Number of genes PathFX identified
    v1genes = data_frame_v1.Gene.shape[0]
    v2genes = data_frame_v2.Gene.shape[0]

    if (v1genes > v2genes) & print_stuff:
        print("PathFX v1 identified " + str(v1genes - v2genes) + " more genes than PathFX v2.")
    else:
        if (v2genes > v1genes) & print_stuff:
            print("PathFX v2 identified " + str(v2genes - v1genes) + " more genes than PathFX v1.")
        else:
            if print_stuff:
                print("Both versions have " + str(v1genes) + " genes.")

    return v2genes - v1genes


def read_paths():
    # Check if the configuration file exists
    locs = [""] * 2
    temp_locs = locs.copy()

    try:
        with open('pathfx_locs.txt', 'r+') as pathfx_locs:
            locs = pathfx_locs.readlines()

            # If the file is created, but empty
            if len(locs) == 0:
                locs = [""] * 2

            # There is only one line in the config file
            if len(locs) == 1:
                locs = [locs[0] + "\n", ""]

            # Throw away all other data besides the first two lines
            locs = [locs[0], locs[1]]

            temp_locs = locs.copy()

            # Check to make sure dirs contain the script file

            path_incorrect = True

            while path_incorrect:
                try:
                    if 'phenotype_enrichment_pathway.py' in os.listdir(locs[0].strip()):
                        path_incorrect = False

                    # Input is a path but it doesn't contain the file
                    else:
                        print("The path doesn't have any version of PathFX, try again for v1: \n")
                        locv1 = input("Enter the path to PathFX version 1: \n")
                        print("\n")
                        locv1 = locv1.strip()
                        locv1 = locv1 + "/scripts"
                        locs[0] = locv1 + "\n"

                # Input isn't a path
                except IOError:
                    locv1 = input("The config doesn't have a valid path for PathFX v1, enter a path: \n")
                    print("\n")
                    locv1 = locv1.strip()
                    locv1 = locv1 + "/scripts"
                    locs[0] = locv1 + "\n"

            path_incorrect = True

            while path_incorrect:
                try:
                    if 'phenotype_enrichment_pathway.py' in os.listdir(locs[1].strip()):
                        path_incorrect = False

                    # Input is a path but it doesn't contain the file
                    else:
                        locv2 = input("The path doesn't have any version of PathFX, try again for v2: \n")
                        print("\n")
                        locv2 = locv2.strip()
                        locv2 = locv2 + "/scripts"
                        locs[1] = locv2

                # Input isn't a path
                except IOError:
                    locv2 = input("The config doesn't have a valid path for PathFX v2, enter a path: \n")
                    print("\n")
                    locv2 = locv2.strip()
                    locv2 = locv2 + "/scripts"
                    locs[1] = locv2

    except IOError:
        # File doesn't exist, so create and then run this function again
        f = open('pathfx_locs.txt', 'w')
        f.close()
        read_paths()

    if locs != temp_locs:
        # Delete current contents and write the new paths
        f = open('pathfx_locs.txt', 'w')
        f.write(locs[0].strip() + "\n" + locs[1])

    return [locs[0].strip(), locs[1].strip()]


def compare_n_drugs(num_of_drugs):
    paths = read_paths()

    # Path to database of all drugs
    drug_db_path = paths[1][0:len(paths) - 10] + "/rscs/repoDB_082117.txt"

    # Read into data frame with pandas
    drug_db = pd.read_csv(drug_db_path)

    # Select n drugs randomly
    random_drugs = drug_db.sample(n=num_of_drugs)

    # Init array with space for n drugs
    selected_drug_names = [""] * num_of_drugs
    # Counter
    i = 0

    # Save original working dir before changing it
    ogwd = os.getcwd()

    # Print progress to the user
    prog = 1

    # Run analyses on the randomly selected drugs
    for drug in random_drugs.drug_id:
        # Change working dir so the script can find its resources
        os.chdir(paths[0])
        os.system("python3 phenotype_enrichment_pathway.py -d " + drug +
                  " -a PathFX_analysis")

        os.chdir(paths[1])
        os.system("python3 phenotype_enrichment_pathway.py -d " + drug +
                  " -a PathFX_analysis")

        print("\n\n" + str(prog) + " out of " + str(num_of_drugs) + "\n")
        prog += 1

    # Set working dir back to original dir
    os.chdir(ogwd)

    # Keep track of how many more genes PathFX discovered.
    # Number will be negative if v1 identified more genes
    gene_diff = 0

    for drug in random_drugs.drug_id:
        v1db = paths[0][0:len(paths) - 10] + "/results/pathfx_analysis/" + drug + "/" + drug + \
               "_merged_neighborhood__assoc_database_sources_.txt"

        v2db = paths[1][0:len(paths) - 10] + "/results/pathfx_analysis/" + drug + "/" + drug + \
            "_merged_neighborhood__assoc_database_sources_.txt"

        gene_diff += compare_one_drug(v1db, v2db, print_stuff=False)

    # Calculate average of how many more genes v2 identified
    gene_diff /= num_of_drugs

    # Compile array for drug name and drug bank ID's
    drug_names = [""] * num_of_drugs
    count = 0

    for drug in random_drugs.drug_name:
        drug_names[count] = drug
        count += 1

    count = 0

    for drug in random_drugs.drug_id:
        drug_names[count] = drug_names[count] + " (" + drug + ")"
        count += 1

    print("\nPathFX v2 identified " + str(gene_diff) + " more genes on average for drugs: " + str(drug_names) + "\n")

    return gene_diff


if args.random_drugs:
    compare_n_drugs(args.random_drugs)

# if len(sys.argv) == 1:
#    compare_one_drug("", "")

# if len(sys.argv) == 2:
#    try:
#        # Make sure arg is a number
#        compare_n_genes(int(sys.argv[1]))

#    except ValueError:
#        print("Please enter a number as the argument.\n")

# if len(sys.argv) == 3:
#    compare_one_drug(sys.argv[1], sys.argv[2])

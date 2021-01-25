# Author: Kyle Petkovic
from functools import partial
import pandas as pd
import sys
import os
import argparse
from multiprocessing import Pool

# Dictates how many data points to show for the most common points when only doing one drug
top_data = 10

print("\n")

# Setting up argparse
parser = argparse.ArgumentParser(description="Calculate differences between between identified phenotypes and genes "
                                             "to help compare PathFX versions")

command_group_optional_args = parser.add_mutually_exclusive_group()
command_group_optional_args_analysis = parser.add_mutually_exclusive_group()
command_group_optional_args.add_argument('-d', '--drugs', metavar='drug_list', type=str, nargs='+',
                                         help='List of drugs to compare with PathFX versions. Compares '
                                              'genes by default. Use --phenotypes or --genes-with-phenotype to change '
                                              'this.')
command_group_optional_args.add_argument('-r', '--random-drugs', metavar="n_drugs", type=int,
                                         help='Number of drugs to randomly compare. Compares genes by default. Use '
                                              '--phenotypes or --genes-with-phenotype to change this.')
command_group_optional_args_analysis.add_argument('-p', '--phenotypes', help='Compares PathFX using the amount of '
                                                                             'phenotypes it identifies',
                                                  action='store_true')
command_group_optional_args_analysis.add_argument('-gwp', '--genes-with-phenotype',
                                                  help='Compares PathFX using only genes that have an associated '
                                                       'phenotype',
                                                  action='store_true')

parser.add_argument('-c', '--cache', help='Use already completed analysis if a previously analyzed drug is analyzed '
                                          'again', action='store_true')

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


# Returns how many more genes v2 identified than v1
def compare_one_drug(v1db, v2db, drug_list, print_stuff=True):
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

    # init values
    v1_data = -1
    v2_data = -1

    # genes with associated phenotypes
    if "_merged_neighborhood__assoc_database_sources_.txt" in v1db:
        data_frame_v1 = pd.read_csv(v1db, delimiter="\t", index_col=False)

        # v2 algorithm does not give the phenotype column a header, so headers have to be overridden
        data_frame_v2 = pd.read_csv(v2db, delimiter="\t", index_col=False, header=0, names=['Gene', 'CUI', 'Phenotype',
                                                                                            'Source Databases'])
        # Number of genes PathFX identified
        v1_data = data_frame_v1.drop_duplicates(subset=['Gene']).shape[0]
        v2_data = data_frame_v2.drop_duplicates(subset=['Gene']).shape[0]

        if len(drug_list) == 1:
            print("\nv1 top 10 most common genes with associated phenotypes\n")
            print(data_frame_v1['Gene'].value_counts()[:top_data])
            print("\n\nv2 top 10 most common genes with associated phenotypes\n")
            print(data_frame_v2['Gene'].value_counts()[:top_data])

    # Phenotypes
    elif "_merged_neighborhood__assoc_table_.txt" in v1db:
        data_frame_v1 = pd.read_csv(v1db, delimiter="\t", index_col=False)
        data_frame_v2 = pd.read_csv(v2db, delimiter="\t", index_col=False)

        # Number of genes PathFX identified
        v1_data = data_frame_v1.drop_duplicates(subset=['phenotype']).shape[0]
        v2_data = data_frame_v2.drop_duplicates(subset=['phenotype']).shape[0]

        # Print top 10 most common data points
        if len(drug_list) == 1:
            print("\nv1 top 10 most common phenotypes\n")
            print(data_frame_v1['phenotype'].value_counts()[:top_data])
            print("\n\nv2 top 10 most common phenotypes\n")
            print(data_frame_v2['phenotype'].value_counts()[:top_data])

    # Genes
    elif "_merged_neighborhood_.txt" in v1db:
        data_frame_v1 = pd.read_csv(v1db, delimiter="\t", index_col=False, names=['Gene', 'col2', 'col3'])
        data_frame_v2 = pd.read_csv(v2db, delimiter="\t", index_col=False, names=['Gene', 'col2', 'col3'])

        # Number of genes PathFX identified
        v1_data = data_frame_v1.drop_duplicates(subset=['Gene']).shape[0]
        v2_data = data_frame_v2.drop_duplicates(subset=['Gene']).shape[0]

        # Print top 10 most common data points
        if len(drug_list) == 1:
            print("\nv1 top 10 most common genes\n")
            print(data_frame_v1['Gene'].value_counts()[:top_data])
            print("\n\nv2 top 10 most common genes\n")
            print(data_frame_v2['Gene'].value_counts()[:top_data])

    if (v1_data > v2_data) & print_stuff:
        print("PathFX v1 identified " + str(v1_data - v2_data) + " more genes than PathFX v2.")
    else:
        if (v2_data > v1_data) & print_stuff:
            print("PathFX v2 identified " + str(v2_data - v1_data) + " more genes than PathFX v1.")
        else:
            if print_stuff:
                print("Both versions have " + str(v1_data) + " genes.")

    return v2_data - v1_data


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


# Run analysis with PathFX
def run_analysis_v1(paths_in, drug_in):
    # Change working dir so the script can find its resources
    os.chdir(paths_in[0])

    # Folder already exists
    if args.cache and os.path.isdir(paths_in[0] + "/../results/pathfx_analysis/" + drug_in):
        # check if text files are there
        files = '\t'.join(os.listdir(paths_in[0] + "/../results/pathfx_analysis/"+drug_in))
        if "_merged_neighborhood_.txt" in files and "_merged_neighborhood__assoc_table_.txt" in files and \
                "_merged_neighborhood__assoc_database_sources_.txt" in files:
            return

    os.system("python3 phenotype_enrichment_pathway.py -d " + drug_in +
              " -a PathFX_analysis")
    # print("\n\n" + str(prog) + " out of " + str(num_of_drugs) + "\n")


# Run analysis with PathFX
def run_analysis_v2(paths_in, drug_in):
    os.chdir(paths_in[1])

    # Folder already exists
    if args.cache and os.path.isdir(paths_in[1] + "/../results/pathfx_analysis/" + drug_in):
        # check if text files are there
        files = '\t'.join(os.listdir(paths_in[0] + "/../results/pathfx_analysis/"+drug_in))
        if "_merged_neighborhood_.txt" in files and "_merged_neighborhood__assoc_table_.txt" in files and \
                "_merged_neighborhood__assoc_database_sources_.txt" in files:
            return

    os.system("python3 phenotype_enrichment_pathway.py -d " + drug_in +
              " -a PathFX_analysis")

    # print("\n\n" + str(prog) + " out of " + str(num_of_drugs) + "\n")


# If num_of_drugs is -1, drug_list must have a passed value
def compare_n_drugs(num_of_drugs=-1, drug_list=None, compare_with_phen=False, compare_gwp=False):
    if drug_list is None:
        drug_list = []

    paths = read_paths()

    # Path to database of all drugs
    drug_db_path = paths[1][0:len(paths) - 10] + "/rscs/repoDB_082117.txt"
    # Read into data frame with pandas
    drug_db = pd.read_csv(drug_db_path)

    if len(drug_list) == 0:
        # Select n drugs randomly
        random_drugs = drug_db.sample(n=num_of_drugs)

        drug_list = random_drugs.drug_id

    # Validate passed in drug list
    else:
        for drug in drug_list:
            if (drug.upper() not in drug_db.drug_id.values) and (drug.title() not in drug_db.drug_name.values):
                print(drug + " is not a valid drug.")
                sys.exit(1)
            num_of_drugs = len(drug_list)

    # Save original working dir before changing it
    ogwd = os.getcwd()

    # Run analyses on the randomly selected drugs
    with Pool(20) as pool:
        analysis_v1 = partial(run_analysis_v1, paths)
        analysis_v2 = partial(run_analysis_v2, paths)
        pool.imap_unordered(analysis_v1, drug_list)
        pool.imap_unordered(analysis_v2, drug_list)
        pool.close()
        pool.join()

    # Set working dir back to original dir
    os.chdir(ogwd)

    # Compares the output depending on the commandline argument
    if compare_with_phen:
        gene_diff = compare_with_phenotypes(drug_list, paths)
        output_comparison = " more phenotypes"
    elif compare_gwp:
        gene_diff = compare_with_genes_with_phenotypes(drug_list, paths)
        output_comparison = " more genes with associated phenotypes"
    else:
        gene_diff = compare_with_genes(drug_list, paths)
        output_comparison = " more genes"

    # Compile array for drug name and drug bank ID's
    drug_names = [""] * num_of_drugs
    count = 0

    for drug in drug_list:
        drug_names[count] = drug
        count += 1

    # Just one drug, so say how much more data v2 found
    if len(drug_list) == 1:
        print("\nPathFX v2 identified " + str(gene_diff) + output_comparison + " for the drug: " + str(drug_names) +
              "\n")
        return gene_diff

    # Calculate average of how many more genes v2 identified
    gene_diff /= num_of_drugs

    print("\nPathFX v2 identified " + str(gene_diff) + output_comparison + " on average for the drugs: " +
          str(drug_names) + "\n")

    return gene_diff


def compare_with_genes(drug_list, paths):
    gene_diff = 0

    for drug in drug_list:
        v1db = paths[0][0:len(paths) - 10] + "/results/pathfx_analysis/" + drug + "/" + drug + \
               "_merged_neighborhood_.txt"

        v2db = paths[1][0:len(paths) - 10] + "/results/pathfx_analysis/" + drug + "/" + drug + \
            "_merged_neighborhood_.txt"

        gene_diff += compare_one_drug(v1db, v2db, drug_list, print_stuff=False)

    return gene_diff


def compare_with_phenotypes(drug_list, paths):
    gene_diff = 0

    for drug in drug_list:
        v1db = paths[0][0:len(paths) - 10] + "/results/pathfx_analysis/" + drug + "/" + drug + \
               "_merged_neighborhood__assoc_table_.txt"

        v2db = paths[1][0:len(paths) - 10] + "/results/pathfx_analysis/" + drug + "/" + drug + \
            "_merged_neighborhood__assoc_table_.txt"

        gene_diff += compare_one_drug(v1db, v2db, drug_list, print_stuff=False)

    return gene_diff


def compare_with_genes_with_phenotypes(drug_list, paths):
    # Keep track of how many more genes PathFX discovered.
    # Number will be negative if v1 identified more genes

    gene_diff = 0

    for drug in drug_list:
        v1db = paths[0][0:len(paths) - 10] + "/results/pathfx_analysis/" + drug + "/" + drug + \
               "_merged_neighborhood__assoc_database_sources_.txt"

        v2db = paths[1][0:len(paths) - 10] + "/results/pathfx_analysis/" + drug + "/" + drug + \
            "_merged_neighborhood__assoc_database_sources_.txt"

        gene_diff += compare_one_drug(v1db, v2db, drug_list, print_stuff=False)

    return gene_diff


if args.random_drugs:
    compare_n_drugs(num_of_drugs=args.random_drugs, compare_with_phen=args.phenotypes,
                    compare_gwp=args.genes_with_phenotype)
elif args.drugs:
    compare_n_drugs(drug_list=args.drugs, compare_with_phen=args.phenotypes,
                    compare_gwp=args.genes_with_phenotype)

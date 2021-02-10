# benchmark_PathFX
Benchmarking performance on PathFX versions

A repository for code comparing PathFX outputs between versions, and eventually input data.

# Setup

Uses the same environment as [PathFX](https://github.com/jenwilson521/PathFX).

Create a python 3.6 virtual environment (E.g. `conda create --name pathfxenv python=3.6`), activate it (`conda activate pathfxenv`), and then use `pip3 install -r requirements.txt` to install the required dependencies. Tested on miniconda 4.9.2.

# Usage

Running the script with no arguments or -h will give an overview of each of the arguments. The first time you run the script it will ask for paths to v1 and v2 of [PathFX](https://github.com/jenwilson521/PathFX). To get version 1.0 and 2.0 clone the respective branches.

### Compare specific drugs (--drugs)

Running compare_versions.py with -d or --drugs and then a space seperated list of drugs will compare PathFX versions using the amount of genes each version identifies. It compares genes by default. For multi-word drugs put it in quotes. Examples: `python compare_versions.py -d metformin` and `python compare_versions.py -d metformin 'Quinidine barbiturate'`

- If only one drug is compared, data about the frequency of each gene/phenotype/gene with associated phenotype will be displayed, otherwise the script will find the average of how many more genes/phenotypes/gene with associated phenotypes that PathFX v2 found.

Example output for running `python compare_versions.py -d metformin`:


Example output for running `python compare_versions.py -d metformin 'Quinidine barbiturate'`:

### Compare random drugs (--random-drugs)

To compare random drugs use --random-drugs or -r and then a number. Example usage: `python compare_versions.py -r 2` will compare 2 random drugs.

- If only one drug is compared, data about the frequency of each gene/phenotype/gene with associated phenotype will be displayed, otherwise the script will find the average of how many more genes/phenotypes/gene with associated phenotypes that PathFX v2 found.

### Compare specific phenotypes (--compare-specific-phenotypes)

To compare a single drug while looking at specific phenotypes use --compare-specific-phenotypes or -csp. Needs to also be used with --drugs. Example usage: `python compare_versions.py -d metformin -csp "lipoprotein(a)"`. Note that phenotype names need to be in quotes if they contain any special characters like spaces or parentheses. 

Example with `python compare_versions.py -d metformin -csp "lipoprotein(a)"`: 

### Compare different data

- genes are compared by default
- to compare phenotypes use --phenotypes or -p
- to compare genes with associated phenotypes use --genes-with-phenotype or -gwp

### Used cached analyses (--cache)

If you know you have already analyzed a drug you can use --cache or -c to use the already computed data. Example: `python compare_versions.py -d metformin -c`. You can also use this when comparing random drugs, but it's unlikely you'll get the same drug over multiple random selections.

### Top amount (--top-amount)

When comparing a single drug by default only the top 10 most frequent genes/phenotypes/genes with associated phenotypes will be displayed. To change this you can run --top-amount or -t and then the amount of data you want to display. Example: `python compare_versions.py -d metformin -t 20` will display the top 20 most frequent genes.


# DODGE: Dynamic Outbreak Detection for Genomic Epidemiology

DODGE is an algorithm and pipeline that identifies potential point source outbreak clusters in bacterial pathogens (denoted investigation clusters) from large scale ongoing genomic surveillance datasets (Allele profiles from cgMLST or SNP calls). Initial clusters should be defined from a background dataset that should ideally represent existing clusters in the population being surveilled. These initial clusters are used as input into the cluster detection script proper.

## Installation

### with conda/mamba

`conda install -c bioconda dodge`

### without conda

`git clone https://github.com/LanLab/dodge.git`

`cd dodge`

`python setup.py install`

**Note following dependencies must be installed**

* scikit-learn
* pandas

## Inputs

### Genetic difference data (cgMLST allele profiles OR SNP data from snippy)

`-i` OR `--variant_data`

**Allele profiles:**  cgMLST allele profiles downloaded from either Enterobase or MGTdb

**SNP calls:** A folder containing 2 files produced by snippy for each isolate in the analysis (strainN.subs.vcf and strainN.consensus.subs.fa) 



### Strain metadata

`-s` OR `--strainmetadata`

**MGTdb:** Strain metadata from MGTdb can be used directly as input into DODGE

**Enterobase:** Strain metadata from Enterobase with hierCC experimental metadata included can be used directly as input into DODGE

**SNP:** For SNP data users must generate a tab delimiteduser metadata file with the following columns:

- 'Strain' or 'Isolate' (for strain identifier)
- 'Year'
- 'Month'
- 'Date' (only required if analysis is weekly)

### Optional distance matrix

`-d` OR `--distances`

To improve performance when analysing large datasets pairwise distance matrices can be precomputed with the pairwise_dist.py script and used as an additional input. Pairwise distances for any isolate with genetic difference data but not present in the distance matrix will be computed.

Tab delimited pairwise distance matrix format :

            StrainA	StrainB	StrainC
    StrainA 0	7	2
    StrainB	7	0	3
    StrainC 2	3	0

### Optional previous clusters

`-c` OR `--inclusters`

A clusters file will be produced by either a --background_data run or a normal run on a given time segment. This cluster file must be used as input for the subsequent time period. i.e. the cluster file from a background data run would be used as input into the first time period of a normal run and the cluster output of one time period run would be used as input into the subsequent time period

## Outputs

### _investigation_clusters.txt file

Tab delimited file containing all information on assigned investigation clusters

**columns as follows:**
* ID - internal arbitrary id used within DODGE
* mgtid	- ID reported for outbreak clusters using nomenclature (for MGT or hierCC)
* Level	- genetic threshold of cluster
* Size - number of isolates assigned
* Max distance - maximum pairwise distance of two isolates within the cluster
* Timespan - number of days or months that the isolates in the cluster span
* Mindate - earliest isolation date of an isolate within the cluster
* Maxdate - latest isolation date of an isolate within the cluster	
* Strains - comma separated list of isolate names
* status - status of that cluster (new, expanded, unchanged)
* Investigation	- Boolean. Whether cluster has been identified as investigation


### _all_clusters.txt file

**Tab delimited file containing all necessary information from existing clusters (investigation and not) from all isolates in the current run**

**columns are identical to _investigation_clusters.txt with the addition of:**
* contains - clusters with lower genetic threshold level that this cluster contains
* partof - clusters with higher genetic threshold level that this cluster is a part of

### _isolate_information.txt file

**Tab delimited file containing information for each isolate.**
**columns are identical to input metadata file with the addition of the following:**

* 0cluster - internal cluster id that this strain is in at 0 genetic distance threshold
* 1cluster - internal cluster id that this strain is in at 1 genetic distance threshold
* 2cluster - internal cluster id that this strain is in at 2 genetic distance threshold
* ...
* Ncluster - internal cluster id that this strain is in at Nth genetic distance threshold
* investigation cluster - When an isolate is in an investigation cluster the ID using nomenclature (for MGT or hierCC)

### _pairwise_distances.txt file

Tab delimited pairwise distance file for all isolates included in the current run matrix format :

            StrainA	StrainB	StrainC
    StrainA 0	7	2
    StrainB	7	0	3
    StrainC 2	3	0



## Full usage:

`dodge.py -i VARIANT_DATA --inputtype {snp,allele} -s STRAINMETADATA --outputPrefix OUTPUTPREFIX [...]`

`-h`, `--help `           show help

**Required input/output:**

`-i`,`--variant_data`: file containing allele profiles (tab delimited table) or snp data (wildcard path to snippy outputs e.g. /folder/*_snippy = /folder/straina_snippy +
                        /folder/strainb_snippy + ...)If using wildcards in path make sure to add "" (default: None)

`--inputtype`: is input data alleles or snps (snp or allele)

`-s`, `--strainmetadata`: file containing isolate information (downloaded from mgtdb, Enterobase or created for SNPs)

`--outputPrefix` output path and prefix for output file generation

**Optional input/output:**

`-d`,` --distances` file containing pairwise distances corresponding to the alleleprofiles file (from previous run of this script if applicable)
  
`-c`,`--inclusters` existing clusters to be imported

  `--background_data` data in this input set / time window to be used for background (no outbreak predictions) (default: False)

  `-n`, `--no_cores` number cores to increase pairwise distance speed (default: 8)

**SNP input specific:**

  `--useref`              include reference in distances/clusters for snp inputtype (default: False)

  `--mask MASK`           bed file for reference used to generate SNPs with regions to ignore SNPs (i.e. phages etc) (default: None)

  `--snpqual SNPQUAL`     minimum allowable SNP quality score (default: 1000)

**Allele input specific:**

  `--enterobase_data`     metadata and allele profiles downloaded from enterobase, if hierCC in metadata table hierCC will be used for outbreak naming (i.e. column named HCXXX)
                        (default: False)

**Date / time options:**

  `--startdate` start date for new cluster analysis (format YYYY-MM-DD if timesegment = week or YYYY-MM if timesegment = month) if left blank earliest date not in inclusters will be identified from strain metadata 

  `--enddate` end date for new cluster analysis (format YYYY-MM-DD) if left blank latest date in input metadata will be used (default: None)
  
  `--timesegment` time segment to perform analysis. every month or every week (default: week)
  
`-t`, `--timewindow` time period a cluster must fall into to be called as investigation --outbreakmethod dodge only --timesegment week default 28
                        --timesegment month default 2 (default: None)

**Clustering options:**

  `-l`, `--dist_limits` comma separated list of cluster cutoffs or range or both i.e `1,2,5` or `1-8` or `1,2,5-10` (default: 1-5)
  
`-m`, `--max_missmatch` maximum number of missmatches reported between 2 isolates (will default to max of --dist_limits + 1 if not set) (default: 1)

**Outbreak detection algorithm options:**

  `--minsize` smallest cluster size for outbreak detection (default: 5)

  `--outbreakmethod` algorithm for outbreak detection dodge or static (default: dodge)

  `--static_cutoff` cutoff for static genetic cutoff method, must be used with `--outbreakmethod static` (default: 5)
  
## Examples

### Example 1 - Retrospective analysis of two months with all available background data (Australian dataset from paper)

For this example all data is already available in one set of input files (i.e. allele profiles and metadata are present for all isolates)

**Step 1** - run dodge with `--background_data` flag on the 1 year of background data

    dodge -a examples/example_aus_2month_w_background_alleleprofile.txt --inputtype allele -s examples/example_aus_2month_w_background_metadata.txt --outputPrefix /path/to/output_folder/prefix_ -n 8  --enddate 2016-12-31 --timesegment week -t 28 -l 1-5 --outbreakmethod dodge --background_data

This command will identify all isolates befire date specified (`--enddate 2021-12-31`) and produce three outputs:
* prefix_background_pairwise_distances.txt
* prefix_background_all_clusters.txt
* prefix_background_isolate_information.txt

**Step 2** - run dodge on the two months to be analysed

    dodge -a examples/example_aus_2month_w_background_alleleprofile.txt --inputtype allele -s examples/example_aus_2month_w_background_metadata.txt --outputPrefix /path/to/output_folder/prefix -n 8  --startdate 2017-01-01 --enddate 2017-02-28 --timesegment week -t 28 -l 1-5 --outbreakmethod dodge -d background_pairwise_distances.txt -c background_all_clusters.txt

This command will check the number of time segments (`--timesegment week`) that are between the specified start and end dates. For this example 9 weeks at least partially fall within the 2 months specified (`--startdate 2022-01-01` `--enddate 2022-02-28`). dodge will internally run 9 sequential runs on one week of isolates at a time. This will produce 9 sets of 4 files with each week producing:

* _pairwise_distances.txt
* _all_clusters.txt
* _isolate_information.txt
* _investigation_clusters.txt (if any are called)

The final outputs will be the files named with the last of the 9 weeks:
* prefix_2017-02-26_2017-03-04_pairwise_distances.txt
* prefix_2017-02-26_2017-03-04_all_clusters.txt
* prefix_2017-02-26_2017-03-04_isolate_information.txt
* **prefix_2017-02-26_2017-03-04_investigation_clusters.txt**

For most cases the investigation_clusters file from the final week (bold above) will provide sufficient information.

### Example 2 - Prospective analysis of 2025 months with 2023,2024 as background using months

For this example the background dataset (2023 and 2024) is already available in one set of input files (i.e. allele profiles and metadata are present for all isolates) however each month of data in 2025 will be run as it is produced using the previous months output files as inputs.

**Step 1** - run dodge with `--background_data` flag on the 1 year of background data

    dodge -a 2023-2024_allele_profile_file.txt --inputtype allele -s 2023-2024_strain_metadata_file.txt --outputPrefix /path/to/output_folder/prefix -n 8  --startdate 2023-01 --enddate 2024-12 --timesegment month -t 2 -l 1-5 --outbreakmethod dodge --background_data

This command will identify all isolates within the date range specified (`--startdate 2023-01` `--enddate 2024-12`) and produce three outputs:
* background_pairwise_distances.txt
* background_all_clusters.txt
* background_isolate_information.txt

**Step 2** - run dodge on the first month to be analysed with allele profile and metadata files containing the background AND new isolates for the current month

    dodge -a 2023-2024_w_2025-01_allele_profile_file.txt --inputtype allele -s 2023-2024_w_2025-01_strain_metadata_file.txt --outputPrefix /path/to/output_folder/prefix -n 8  --startdate 2025-01 --enddate 2025-01 --timesegment month -t 2 -l 1-5 --outbreakmethod dodge -d background_pairwise_distances.txt -c background_all_clusters.txt

This command will identify all isolates within the date range specified (`--startdate 2023-01` `--enddate 2024-12`) and produce three outputs:
* prefix_2025-01_pairwise_distances.txt
* prefix_2025-01_all_clusters.txt
* prefix_2025-01_isolate_information.txt

If investigation clusters are identified then a fourth file will be generates

* prefix_2025-01_investigation_clusters.txt

**Steps 3 to 25** - run dodge on the subsequent months to be analysed with allele profile and metadata files containing the background, previous months AND new isolates for the current month

    dodge -a background+current_and_previous_months_allele_profile_file.txt --inputtype allele -s background+current_and_previous_months_strain_metadata_file.txt --outputPrefix /path/to/output_folder/prefix -n 8  --startdate 2025-02 --enddate 2025-02 --timesegment month -t 2 -l 1-5 --outbreakmethod dodge -d prefix_2025-01_pairwise_distances.txt -c prefix_2025-01_all_clusters.txt

Outputs will be as per the first month but named for the month they represent.

## Generating pairwise distances only with dodgedists

This script only runs the distance matrix generation component of dodge. Inputs files are the same as dodge except no strain metadata file is required. Output of this script can be passed to dodge to reduce running time of the main dodge script.

`dodgedists [-h] -i VARIANT_DATA --inputtype {snp,allele} --output OUTPUT [-d DISTANCES] [-n NO_CORES] [-m MAX_MISSMATCH] [--useref] [--mask MASK] [--snpqual SNPQUAL]
                        [--enterobase_data]`

`-h`, `--help `           show help

**Required input/output:**

`-i`,`--variant_data`: file containing allele profiles (tab delimited table) or snp data (wildcard path to snippy outputs e.g. /folder/*_snippy = /folder/straina_snippy +
                        /folder/strainb_snippy + ...)If using wildcards in path make sure to add "" (default: None)

`--inputtype`: is input data alleles or snps (snp or allele)

`--output` output path and prefix for output file generation

**Optional input/output:**

`-d`,` --distances` file containing pairwise distances corresponding to the alleleprofiles file (from previous run of this script if applicable)
  
**Run options:**

`-n`, `--no_cores` number cores to increase pairwise distance speed (default: 8)

`-m`, `--max_missmatch` maximum number of missmatches reported between 2 isolates (will default to max of --dist_limits + 1 if not set) (default: 1)

**SNP input specific:**

  `--useref`              include reference in distances/clusters for snp inputtype (default: False)

  `--mask MASK`           bed file for reference used to generate SNPs with regions to ignore SNPs (i.e. phages etc) (default: None)

  `--snpqual SNPQUAL`     minimum allowable SNP quality score (default: 1000)

**Allele input specific:**

  `--enterobase_data`     metadata and allele profiles downloaded from enterobase, if hierCC in metadata table hierCC will be used for outbreak naming (i.e. column named HCXXX)
                        (default: False)

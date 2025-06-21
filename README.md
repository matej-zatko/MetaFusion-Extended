# MetaFusion-Extended

MetaFusion-Extended is an extension of [MetaFusion Clinical](https://github.com/ccmbioinfo/MetaFusion-Clinical) that adds enhanced gene symbol unification and graph-based clustering using gene identifiers to improve the merging of fusion calls across multiple tools. 

In both the original MetaFusion and MetaFusion Clinical, gene renaming relies on the NCBI database to update gene names to their official symbols. However, this step can fail if a gene symbol is missing from the database, leading to unsuccessful merging based on gene names. To mitigate this problem, MetaFusion-Extended integrates multiple gene symbol databases: [NCBI](https://www.ncbi.nlm.nih.gov/gene/), [HGNC](https://www.genenames.org/), and [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html). Additionaly, gene symbols are not stable and can change or even be removed over time. Therefore, MetaFusion-Extended also incorporates Ensembl gene IDs into its graph clustering step (see the [original MetaFusion manuscript](https://academic.oup.com/bioinformatics/article/37/19/3144/6263829)) to improve fusion call merging across different tools. 

Note that one limitation of this graph-based clustering approach is the inability to detect fusion isoforms. If you are interested in isoform detection, [MetaFusion Clinical](https://github.com/ccmbioinfo/MetaFusion-Clinical), which only uses clusters by breakpoints, is recommended instead.

## Installation

The recommended way to run MetaFusion-Extended is via Docker, using a pre-built image with all dependencies included. For alternative installation methods (e.g., local installation or Singularity), please refer to the [MetaFusion](https://github.com/ccmbioinfo/MetaFusion/wiki) and [MetaFusion Clinical](https://github.com/ccmbioinfo/MetaFusion-Clinical/wiki) wikis.

### Step 1: Clone the repository

The code can be downloaded by cloning the repository using the following command:

```
git clone https://github.com/matej-zatko/MetaFusion-Extended.git
```

This will create a folder named MetaFusion-Extended in the directory where the above command was executed.

### Step 2: Pull the Docker image

Download the Docker image with all dependencies pre-installed:

```
docker pull mapostolides/metafusion:readxl_writexl
```

### Step 3: Download reference files

Reference files have to be downloaded separately from [here](https://drive.google.com/file/d/1pxKYmG3LYOccdJWnfhDZ-LDARuMiace6/view?usp=sharing), unzipped and placed inside the MetaFusion-Extended folder (`MetaFusion-Extended/reference_files`).

### Step 4: Run the container

Start an interactive Docker container with your project directory mounted (replace `/your/local/path/MetaFusion-Extended` with the absolute path to your cloned repository):

```
docker run --name MetaFusion-Extended --rm -it --entrypoint /bin/bash -v /your/local/path/MetaFusion-Extended:/MetaFusion mapostolides/metafusion:readxl_writexl
```

You will now be inside the container and ready to run MetaFusion-Extended commands from within the `/MetaFusion` directory.

## Usage

### Generating CFF files

Before running MetaFusion-Extended, the output from individual fusion callers must be converted to CFF (Common Fusion Format) and concatenated into a single file. MetaFusion-Extended extends this format by adding fields containing Ensembl IDs. To learn more about the format, see [Metafusion file formats](https://github.com/ccmbioinfo/MetaFusion/wiki/metafusion-file-formats).

To generate the CFF file, use the following script:

```
/MetaFusion/cff_generation/RUN_convert_fusion_results_to_cff.sh \
  -c <caller_file_dir> \
  -s <sampleinfo> \
  -d <dataset> \
  -o <outdir> \
  -t "<tool1 tool2 ...>"
```

Parameters:
* `-c`: Path to the directory containing output files from the fusion callers. This directory must follow an exact structure (see the [MetaFusion wiki](https://github.com/ccmbioinfo/MetaFusion/wiki/How-to-generate-a-CFF-file)). An example of correct directory structure can be seen in `test_data/caller_output_files`.
* `-s`: Path to a sample_info.txt file containing sample metadata. An example can be found in `test_data/caller_output_files/sample_info.txt`
* `-d`: Name of the dataset as it appears in the directory strcuture.
* `-o`: Output directory where the combined CFF file will be written.
* `-t`: Space-separated list of tools used (e.g., `"arriba star_fusion"`). The supported tools are: arriba, star_fusion, star_seqr, defuse, ericscript, integrate, fusionmap, cicero, fusioncatcher.

This will convert the outputs of each tool for each sample and merge them into a single `merged.cff` file, that is then ready to be used as input to MetaFusion-Extended.

For all benchmarking datasets included in `test_data` (see [MetaFusion wiki](https://github.com/ccmbioinfo/MetaFusion/wiki/Benchmarking-dataset-fastq-files) for information about the datasets), extended CFF files have been already generated and are available in the `test_data/cff` folder.

#### Example of generating a CFF file for the BRCA dataset:

```
/MetaFusion/cff_generation/RUN_convert_fusion_results_to_cff.sh -c /MetaFusion/test_data/caller_output_files -s /MetaFusion/test_data/caller_output_files/sample_info.txt -d BRCA -o /MetaFusion/test_data/cff/BRCA -t "arriba defuse ericscript fusionmap integrate star_fusion star_seqr"
```

### Running MetaFusion-Extended

To run MetaFusion-Extended. navigate to the `scripts` directory (important!) and run:

```
bash MetaFusion.extended.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --annotate_exons \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --gene_info $gene_info \
                 --hgnc_db $hgnc_db \
                 --num_tools=2  \
                 --per_sample \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools \
		             --database $database \
                 --update_hist \
                 --ref_dir $ref_dir
```

The input parameters are identical to MetaFusion Clinical and are explained [here](https://github.com/ccmbioinfo/MetaFusion-Clinical/wiki#running-metafusion-clinical). The only additional parameter is `hgnc_db` which specifies the HGNC database file included in `reference_files`. The `per_sample` parameter was removed and fusions in different samples are now on separate lines by default.

This will generate many files containing the intermediate steps of the pipeline. The final output file is `final.n2.cluster.filt` in text format, or `final.n2.cluster.filt.xlsx` as an Excel table. Description of individual columns in the final cluster file is available [here](https://github.com/ccmbioinfo/MetaFusion-Clinical/wiki#final-output-file-column-descriptions).

#### Example of running MetaFusion-Extended with the BRCA CFF generated earlier:

```
bash MetaFusion.extended.sh --outdir /MetaFusion/outputs/BRCA \
                 --cff /MetaFusion/test_data/cff/BRCA.cff  \
                 --gene_bed /MetaFusion/reference_files/new_bed.total.Oct-1-2020.uniq.bed \
                 --annotate_exons \
                 --fusion_annotator \
                 --genome_fasta /MetaFusion/reference_files/human_g1k_v37_decoy.fasta \
                 --gene_info /MetaFusion/reference_files/Homo_sapiens.gene_info.new \
                 --hgnc_db /MetaFusion/reference_files/hgnc_complete_set_2025-01-06.txt \
                 --num_tools=2  \
                 --recurrent_bedpe /MetaFusion/reference_files/blocklist_breakpoints.bedpe \
                 --scripts /MetaFusion/scripts \
                 --database /MetaFusion/reference_files/historical_database.db \
                 --update_hist \
                 --ref_dir /MetaFusion/reference_files
```

### Benchmarking

The performance of MetaFusion-Extended can be benchmarked using the included FusionBenchmarking toolkit. By comparing the final output with a truth set of fusions known to be in the sample, the toolkit will mark fusion calls as either true positives (TP), false positives (FP) or false negatives (FN). More information about the toolkit can be found [here](https://github.com/ccmbioinfo/MetaFusion/wiki/benchmarking_toolkit).

Run the benchmarking toolkit using the following command:

```
bash benchmarking_cluster.MetaFusion.sh <out_dir> <truth_set> <final_cluster> <scripts_dir>
```

Arguments:
* `<out_dir>`: Path to the directory where output files will be written.
* `<truth_set>`: Path to the truth set file.
* `<final_cluster>`: Path to the final cluster file (e.g., the MetaFusion output).
* `<scripts_dir>`: Path to the scripts directory.

Original truth sets for the included datasets are available in the `test_data/truth_sets` directory. Updated truth sets, in which gene symbols have been updated using the enhanced gene renaming script (`scripts/update_truth_set_gene_symbols.R`), can be found in `test_data/updated_truth_sets`.

#### Example of running the benchmarking toolkit on the BRCA cluster file:

```
bash benchmarking_cluster.MetaFusion.sh /MetaFusion/outputs/BRCA/FusionBenchmarking /MetaFusion/test_data/updated_truth_sets/BRCA.truth_set.dat /MetaFusion/outputs/BRCA/final.n2.cluster.filt /MetaFusion/scripts
```

## Links to original implementations

### MetaFusion-Clinical

MetaFusion Clinical implementation can be found here: https://github.com/ccmbioinfo/MetaFusion-Clinical. 

The associated manuscript was published in the Journal of Molecular Diagnostics: https://www.jmdjournal.org/article/S1525-1578(23)00211-8/fulltext.

### Initial implementation of MetaFusion

Initial implementation of MetaFusion can be found here: https://github.com/ccmbioinfo/MetaFusion. 

The manuscript published in Bioinformatics: https://academic.oup.com/bioinformatics/article/37/19/3144/6263829

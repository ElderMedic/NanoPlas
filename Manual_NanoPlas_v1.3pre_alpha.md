# Nanopore Plasmid Sequencing Pipeline

This pipeline is designed for processing Nanopore sequencing data from plasmids. It performs various steps, including basecalling, demultiplexing, quality control, assembly, annotation, variant calling, and result generation.

## Usage

To run the pipeline, execute the following command:

```bash
NanoPlas --working-dir <working_directory> [--config <config_file>] [--sample-config <sample_config_file>] [--no-guppy] [--skip <steps_to_skip>] [--quiet]
```

Location of the pipeline scripts: `/home/macrogen-ont/NanoPlas/`

Replace `<working_directory>` with the path to the base directory where the output files will be generated. Optionally, you can provide the names (not path) of the configuration file (`<config_file>`) and sample configuration file (`<sample_config_file>`), but they have to present at the `<working_directory>`. By default, the pipeline looks for `pipeline_config.json` and `sample_config.tsv` in the `<working_directory>`.

To fill in the best-practice configuration files, please refer to the [Standard config](./standard_config/) folder.

### Optional Arguments

- `--no-guppy`: Skip the Guppy basecalling and barcoding steps.
- `--skip <steps_to_skip>`: Skip the specified pipeline steps. Provide a space-separated list of steps to skip. Available steps are: `mergefq`, `assembly`, `plannotate`, `minimap2`, `samtools`, `bcftools`, `variant_call`, `copy`.
- `--quiet`: Run the pipeline in quiet mode, suppressing stdout/stderr debug info of cmd.

## Dependencies

Make sure you have the following dependencies installed before running the pipeline:

conda env: base (temp), will be moved to another env.

- Guppy
- Flye
- Medaka
- Sniffles
- Samtools
- Bcftools
- Plannotate
- NanoPlot
- Minimap2
- Tabix

Refer to the documentation of each tool for installation instructions and additional requirements.

## Configuration

The pipeline requires a configuration file (`pipeline_config.json`) and a sample configuration file (`sample_config.tsv`). These files contain the necessary settings and information for running the pipeline.

The `pipeline_config.json` file includes various parameters for each step of the pipeline, such as input/output directories, quality thresholds, and tool-specific settings. Settings in this config is consistent and fixed across samples.

The `sample_config.tsv` file provides information about each sample, including the assigned barcodes, plasmid sizes, and other sample-specific details. This will be used during iteration per sample.

Make sure to customize the configuration files according to your specific experiment and requirements before running the pipeline and after for potential reruns.

## Workflow

The pipeline follows the following workflow: all steps are optional but you need to include at least one of them.

1. Guppy basecalling and demultiplexing .
2. Quality control of basecalled fastq files, visualize and generate report with Nanoplot.
3. Assembly using Flye, polish with medaka.
4. Annotation using Plannotate.
5. Alignment using Minimap2.
6. BAM file processing using Samtools: sort, rmdup, index.
7. Variant calling using Bcftools. _bcftools.vcf.gz is final indexed 
8. Variant calling using Medaka haploit variant.
9. Structural variant calling using Sniffles2.
10. Convert vcf files to csv format.
11. Copying result files to the output directory.

The pipeline performs these steps for each barcode present in the `sample_config.tsv` file. Intermediate output and results are saved in the specified output directory (`<working_directory>/Barcoding`).

## Output

The pipeline generates various output files and directories, including:

- Basecalled and demultiplexed fastq files, merged fastq file.
- Quality control reports (NanoPlot).
- Assembly results (FASTA files). if applicable
- Annotation results (HTML files). if applicable
- Alignment results (SAM/BAM files). if applicable
- Variant calling results (VCF files). if applicable
  - medaka variant calling results. if applicable
  - bcftools variant calling results. if applicable
  - Sniffles2 Structural variant calling results. if applicable.
  - variant summary statistics csv files. if applicable

The pipeline also copies selected result files to the output directory specified in the configuration.

## Contact

If you have any questions or need further assistance, please contact Changlin or Jongbum.

---

## References

If you use this pipeline in your research, please consider citing the relevant tools and references:

- Guppy: [Oxford Nanopore Technologies](https://nanoporetech.com)
- Flye: [Kolmogorov et al., 2019](https://www.nature.com/articles/s41592-019-0669-3)
- Medaka: [GitHub Repository](https://github.com/nanoporetech/medaka)
- Sniffles: [Sedlazeck et al., 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1611-3)
- Samtools: [GitHub Repository](https://github.com/samtools/samtools)
- Bcftools: [GitHub Repository](https://github.com/samtools/bcftools)
- Plannotate: [GitHub Repository](https://github.com/sanger-pathogens/plannotate)
- NanoPlot: [GitHub Repository](https://github.com/wdecoster/NanoPlot)
- Minimap2: [GitHub Repository](https://github.com/lh3/minimap2)
- Tabix: [GitHub Repository](https://github.com/samtools/tabix)

Please refer to the documentation and publications of these tools for more information.

---

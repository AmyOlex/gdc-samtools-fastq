# gdc-samtools-fastq

A WDL pipeline to convert GDC BAM files to fastq format utilizing GDC recommended options. For execution on Terra.bio.

# GDC Samtools FASTQ Workflow
This WDL workflow converts aligned BAM files back into paired-end FASTQ format, adhering to GDC Data Harmonization standards. It is designed to run locally using the Dockstore CLI or in the cloud on Terra/AnVIL.

## Version History

- **v2.0 (HTTPS Download Edition):** Downloads BAM files directly from the GDC HTTPS API, bypassing DRS/GCS resolution. Use this version if DRS resolution is failing on Terra due to broken GDC GCS buckets. Adds configurable resource settings and preemptible VM support.
- **v1.x (DRS/GCS Edition):** Uses Terra's built-in DRS resolution to access BAM files via Google Cloud Storage. Recommended when GDC GCS buckets are functional.

## Overview

The workflow performs the following steps:
 1) **Downloads BAM (v2.0):** Downloads the BAM file directly from the GDC HTTPS API using your GDC authentication token.
 2) **Indexes BAM:** Automatically generates a .bai index for the input BAM.
 3) **Splits BAM by Read Group:** Detects @RG tags in the BAM header and splits the file into separate BAMs for each read group (lane).
 4) **Restores Original Qualities:** Uses samtools fastq -O to restore original quality scores (OQ tag) if they were recalibrated.
 5) **Parallel Conversion:** Converts each split BAM into paired FASTQ files (R1 and R2) in parallel.
 6) **Merges Output:** Concatenates all split FASTQ files into a single pair (merged_R1.fastq.gz, merged_R2.fastq.gz) for easy downstream use.

## Requirements

 - **WDL Version:** 1.0
 - **Docker Image:** staphb/samtools:1.22 (default)
 - **Executor:** [Cromwell](https://github.com/broadinstitute/cromwell) (via Dockstore CLI or Terra)
 - **GDC Token (v2.0):** A valid authentication token from the [GDC Data Portal](https://portal.gdc.cancer.gov/) (download under your user profile). Tokens expire periodically and must be refreshed.

## Inputs

### v2.0 (HTTPS Download Edition)

Input Name | Type | Description
| :--- | :--- | :--- |
GDC_Samtools_Fastq.sample_id | String | DRS URI (e.g. `drs://dg.4DFC:008b411e-...`) or bare UUID for the GDC BAM file.
GDC_Samtools_Fastq.gdc_token | File | GDC authentication token file for HTTPS downloads.
GDC_Samtools_Fastq.docker_image | String | (Optional) Docker image to use. Default: `staphb/samtools:1.22`.

### v1.x (DRS/GCS Edition)

Input Name | Type | Description
| :--- | :--- | :--- |
GDC_Samtools_Fastq.input_bam | File | The aligned BAM file to convert. Accepts DRS URIs or gs:// paths. Must contain @RG (Read Group) headers.
GDC_Samtools_Fastq.docker_image | String | (Optional) Docker image to use. Default: `staphb/samtools:1.22`.

### Optional Resource Settings (v2.0)

Each task in the workflow has configurable disk, memory, CPU, and preemptible settings. These are all optional and have sensible defaults that work well for typical TCGA BAM files up to ~40 GB. The naming convention is `{task}_{resource}`, for example `download_disk_gb`, `index_memory_gb`, `split_cpu`, or `fastq_preemptible`.

Task | Disk (GB) | Memory (GB) | CPU | Preemptible
| :--- | :---: | :---: | :---: | :---: |
download | 50 | 4 | 2 | 1
index | 50 | 4 | 1 | 3
split | 100 | 4 | 2 | 3
fastq | 50 | 8 | 2 | 3
merge | 100 | 4 | 1 | 3

If your BAM files are larger than ~40 GB, increase the disk sizes accordingly. The split task needs roughly 2x the BAM size, and the merge task needs roughly 2x the total FASTQ output size.

## Outputs
Output Name | Type | Description
| :--- | :--- | :--- |
GDC_Samtools_Fastq.merged_r1 | File | Single merged Read 1 FASTQ file (gzipped).
GDC_Samtools_Fastq.merged_r2 | File | Single merged Read 2 FASTQ file (gzipped).

## Migrating from v1.x to v2.0

1) Upload your GDC token file to your Terra workspace bucket.
2) Set the new `gdc_token` workflow input to point to your uploaded token file.
3) Point `sample_id` to the same data table column you previously used for `input_bam` — DRS URIs work as-is.
4) (Optional) Adjust resource and preemptible settings for your sample sizes. The defaults work well for BAM files up to ~40 GB.

## Understanding Preemptible VMs

Preemptible (spot) VMs are spare Google Cloud machines that cost 60-91% less than standard VMs. The tradeoff is that Google can reclaim them at any time, which kills your running task. Terra will automatically retry.

The preemptible value controls how many times Terra retries on a cheap VM before falling back to a guaranteed standard VM. For example, `preemptible = 3` means try up to 3 times on cheap VMs, then run on a standard VM if all attempts are interrupted.

Guidelines for setting preemptible values:
- Short tasks (under 1 hour): `preemptible = 3` is recommended and saves significant cost.
- The download task: `preemptible = 1` is recommended since a preemption means re-downloading the entire BAM file.
- Long-running tasks (6+ hours): `preemptible = 0` may save time since each retry starts from scratch.

## Running on Terra (v2.0)

 1) Push this repository to GitHub.
 2) Register the workflow on [Dockstore](https://dockstore.org/).
 3) Click "Launch with Terra" if available. If not, manually import using the following URL: [https://app.terra.bio/#import-workflow/dockstore/github.com/AmyOlex/gdc-samtools-fastq/gdc-samtools-fastq:main](https://app.terra.bio/#import-workflow/dockstore/github.com/AmyOlex/gdc-samtools-fastq/gdc-samtools-fastq:main)
 4) Download your GDC token from the [GDC Data Portal](https://portal.gdc.cancer.gov/) and upload it to your workspace bucket.
 5) In your workflow configuration, set `sample_id` to your data table column containing DRS URIs and `gdc_token` to the uploaded token file.
 6) (Optional) Adjust disk, memory, and preemptible settings based on your BAM file sizes.

## Standard Local Testing (Dockstore CLI)

If you have open internet access, use the standard Dockstore method:

1) Install: [Docker Desktop](https://www.docker.com/products/docker-desktop/) and [Dockstore CLI](https://docs.dockstore.org/en/stable/launch-with/launch.html).
2) Run:
   
```bash
dockstore workflow launch \
  --local-entry gdc-samtools-fastq.wdl \
  --json gdc_inputs.json
```
_(Note: For local testing to avoid downloading large files from the web, create a local_inputs.json file and point it to your local BAM file.)_

## Corporate / Restricted Network & Silicon Mac Setup

If you are on a restricted network (VPN) or using an Apple Silicon (M1/M2/M3) Mac where dockstore CLI fails to download files or images, follow this "Offline" procedure, which skips utilizing Dockstore and calls Cromwell directly. This is also a more flexible way to control Cromwell behavior using the cromwell.config file.

### 1. Manual Docker Setup (Silicon Mac)
Pre-pull the image using the specific architecture to ensure compatibility with Terra (Intel/AMD64). This ensures Rosetta handles the translation correctly.

```bash
docker pull --platform linux/amd64 staphb/samtools:1.22
```

### 2. Cromwell Configuration
Create a file named cromwell.conf in your project root. This tells Cromwell not to check Docker Hub for image hashes (which often fails on corporate VPNs).

File: cromwell.conf

```bash
docker {
  hash-lookup {
    enabled = false
  }
}
```

### 3. Generate Local Test Data
Since downloading external BAMs is restricted, run this script to generate a valid, sorted, multi-read-group BAM locally using Docker.

#### Create directory

```bash
mkdir -p test_data
cd test_data
```
#### 1. Create dummy SAM content
```bash
cat <<EOF > test_data.sam
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@RG	ID:Lane1	SM:SampleA	PL:ILLUMINA
@RG	ID:Lane2	SM:SampleA	PL:ILLUMINA
read1_lane1	99	chr1	10	30	10M	=	50	50	AAAAAAAAAA	IIIIIIIIII	RG:Z:Lane1
read1_lane1	147	chr1	50	30	10M	=	10	-50	TTTTTTTTTT	IIIIIIIIII	RG:Z:Lane1
read2_lane2	99	chr1	20	30	10M	=	60	50	GGGGGGGGGG	IIIIIIIIII	RG:Z:Lane2
read2_lane2	147	chr1	60	30	10M	=	20	-50	CCCCCCCCCC	IIIIIIIIII	RG:Z:Lane2
EOF
```

#### 2. Convert to BAM (Sort & Index) via Docker
**Note: We use --platform to match the image we pulled**

```bash
docker run --platform linux/amd64 --rm -v "$PWD":/data -w /data staphb/samtools:1.22 \
  bash -c "samtools view -u test_data.sam | samtools sort -o test_input.bam && samtools index test_input.bam"
```
#### 3. Cleanup text file

```bash
rm test_data.sam
cd ..
```

### 4. Run Directly (Bypassing Dockstore CLI)
Instead of using dockstore workflow launch, run the Cromwell JAR directly. You may need to download the Cromwell JAR manually if the CLI failed to get it.

Run the workflow:

```bash
java -Dconfig.file=cromwell.conf \
  -jar ~/.dockstore/libraries/cromwell-86.jar \
  run gdc-samtools-fastq.wdl \
  --inputs local_inputs.json
```

## Author
Name: Amy Olex

Contact: alolex@vcu.edu

Affiliation: Virginia Commonwealth University

## License
[GNU General Public License v3.0](https://github.com/AmyOlex/gdc-samtools-fastq/blob/main/LICENSE)

## Gen AI Disclaimer
Gemini Pro in Thinking Mode using VCU's secure account was utilized in Dec 2025 to assist in the development of this workflow and writing of documentation. Claude Opus 4.6 via claude.ai was utilized in Feb 2026 to debug the GDC GCS bucket issue and develop the v2.0 HTTPS download edition. All code and documentation was manually reviewed, edited and tested by the author.
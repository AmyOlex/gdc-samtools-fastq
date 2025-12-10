# gdc-samtools-fastq

A WDL pipeline to convert GDC BAM files to fastq format utilizing GDC recommended options.  For execution on Terra.bio.

# GDC Samtools FASTQ Workflow
This WDL workflow converts aligned BAM files back into paired-end FASTQ format, adhering to GDC Data Harmonization standards. It is designed to run locally using the Dockstore CLI or in the cloud on Terra/AnVIL.

## Overview

The workflow performs the following steps:
 - **Splits BAM by Read Group:** Detects @RG tags in the BAM header and splits the file into separate BAMs for each read group (lane).
 - **Restores Original Qualities:** Uses samtools fastq -O to restore original quality scores (OQ tag) if they were recalibrated.
 - **Parallel Conversion:** Converts each split BAM into paired FASTQ files (R1 and R2) in parallel.

## Requirements

 - **WDL Version:** 1.0
 - **Docker Image:** staphb/samtools:1.22 (default)
 - **Executor:** Cromwell (via Dockstore CLI or Terra)

## Inputs

Input Name | Type | Description
| :--- | :--- | :--- |
GDC_Samtools_Fastq.input_bam | File | The aligned BAM file to convert. Must contain @RG (Read Group) headers.
GDC_Samtools_Fastq.input_bam_index | File | The index (.bai) for the input BAM.
GDC_Samtools_Fastq.docker_image | String | (Optional) Docker image to use. Default: staphb/samtools:1.22. Use staphb/samtools:latest for local dev if needed.

## Outputs
Output Name | Type | Description
| :--- | :--- | :--- |
GDC_Samtools_Fastq.r1_fastqs | Array[File] | List of Read 1 FASTQ files (one per Read Group).
GDC_Samtools_Fastq.r2_fastqs | Array[File] | List of Read 2 FASTQ files (one per Read Group).


## Standard Local Testing (Dockstore CLI)

If you have open internet access, use the standard Dockstore method:

1) Install: Docker Desktop and Dockstore CLI.
2) Run:
   
```bash
dockstore workflow launch \
  --local-entry gdc-samtools-fastq.wdl \
  --json gdc_inputs.json
```

## Corporate / Restricted Network & Silicon Mac Setup

If you are on a restricted network (VPN) or using an Apple Silicon (M1/M2/M3) Mac where dockstore CLI fails to download files or images, follow this "Offline" procedure.

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
  --inputs gdc_inputs.json
```

## Running on Terra
 - Push this repository to GitHub.
 - Register the workflow on Dockstore.
 - Click "Launch with Terra".
 - In Terra, upload your BAM files to the workspace bucket and update the inputs to point to gs://... locations.
 - License: MIT
   
```bash
  run gdc-samtools-fastq.wdl \
  --inputs gdc_inputs.json
```

version 1.0

# =============================================================================
# GDC_Samtools_Fastq v2.0 - HTTPS Download Edition
# =============================================================================
# This version downloads BAM files directly from the GDC HTTPS API instead of
# relying on DRS resolution to Google Cloud Storage buckets.
#
# Use this version when:
#   - GDC GCS buckets are unavailable or decommissioned
#   - DRS resolution is failing on Terra
#   - You prefer explicit control over the download method
#
# For the DRS/GCS version (requires working GDC GCS buckets), see v1.x tags.
# =============================================================================


workflow GDC_Samtools_Fastq {
  input {
    # DRS URI or bare UUID for the BAM file
    # Accepts either format:
    #   "drs://dg.4DFC:008b411e-7de8-46f2-9ad8-185cd49ee2e6"
    #   "008b411e-7de8-46f2-9ad8-185cd49ee2e6"
    String sample_id

    # GDC token file for authenticated downloads
    File gdc_token

    # DEFAULT FOR TERRA: Use the public, pinned version
    String docker_image = "staphb/samtools:1.22"

    # Preemptible retry attempts per task (0 = always use standard VMs)
    Int download_preemptible = 1
    Int index_preemptible = 3
    Int split_preemptible = 3
    Int fastq_preemptible = 3
    Int merge_preemptible = 3

    # Resource settings per task
    Int download_disk_gb = 50
    Int download_memory_gb = 4
    Int download_cpu = 2

    Int index_disk_gb = 50
    Int index_memory_gb = 4
    Int index_cpu = 1

    Int split_disk_gb = 100
    Int split_memory_gb = 4
    Int split_cpu = 2

    Int fastq_disk_gb = 50
    Int fastq_memory_gb = 8
    Int fastq_cpu = 2

    Int merge_disk_gb = 100
    Int merge_memory_gb = 4
    Int merge_cpu = 1
  }

  # 0. Download BAM from GDC HTTPS API
  call DownloadFromGDC {
    input:
      sample_id = sample_id,
      gdc_token = gdc_token,
      docker = docker_image,
      disk_gb = download_disk_gb,
      memory_gb = download_memory_gb,
      cpu = download_cpu,
      preemptible = download_preemptible
  }
  # 1. Index BAM file
  call IndexBam {
    input:
      bam = DownloadFromGDC.bam,
      docker = docker_image,
      disk_gb = index_disk_gb,
      memory_gb = index_memory_gb,
      cpu = index_cpu,
      preemptible = index_preemptible
  }
  
  # 2. Split the BAM by Read Group (RG)
  # This is critical for GDC BAMs which often contain multiple lanes
  call SplitBamByRG {
    input:
      bam = DownloadFromGDC.bam,
      bai = IndexBam.bai,
      docker = docker_image,
      disk_gb = split_disk_gb,
      memory_gb = split_memory_gb,
      cpu = split_cpu,
      preemptible = split_preemptible
  }

  # 3. Convert each split BAM into FASTQ pairs (Parallel Scatter)
  scatter (split_bam in SplitBamByRG.split_bams) {
    call BamToFastq {
      input:
        bam = split_bam,
        docker = docker_image,
        disk_gb = fastq_disk_gb,
        memory_gb = fastq_memory_gb,
        cpu = fastq_cpu,
        preemptible = fastq_preemptible
    }
  }

  # 4. Merge the scattered FASTQ files into a single pairs
  call MergeFastqs {
    input:
      r1_files = BamToFastq.r1,
      r2_files = BamToFastq.r2,
      base_name = basename(DownloadFromGDC.bam, ".bam"),
      docker = docker_image,
      disk_gb = merge_disk_gb,
      memory_gb = merge_memory_gb,
      cpu = merge_cpu,
      preemptible = merge_preemptible
  }

  output {
    File merged_r1 = MergeFastqs.merged_r1
    File merged_r2 = MergeFastqs.merged_r2
  }
}


# TASK 0: Download BAM from GDC HTTPS API
task DownloadFromGDC {
  input {
    String sample_id
    File gdc_token
    String docker
    Int disk_gb
    Int memory_gb
    Int cpu
    Int preemptible
  }

  command <<<
    set -e

    # Extract UUID: handle both DRS URIs and bare UUIDs
    UUID=$(echo "~{sample_id}" | sed 's/.*://')
    echo "Resolved UUID: ${UUID}"

    # Download BAM from GDC HTTPS API
    echo "Downloading from https://api.gdc.cancer.gov/data/${UUID}"
    curl -f -L \
      -H "X-Auth-Token: $(cat ~{gdc_token})" \
      -o "${UUID}.bam" \
      "https://api.gdc.cancer.gov/data/${UUID}"

    echo "Download complete. File size:"
    ls -lh "${UUID}.bam"
  >>>

  runtime {
    docker: docker
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible
  }

  output {
    File bam = glob("*.bam")[0]
  }
}


# TASK 1: Index BAM file
task IndexBam {
  input {
    File bam
    String docker
    Int disk_gb
    Int memory_gb
    Int cpu
    Int preemptible
  }

  # Calculate basename so we can symlink with the correct original name
  String filename = basename(bam)

  command <<<
    set -e
    # 1. Symlink input to current working directory
    # This solves the issue of read-only input directories or missing sidecar files
    ln -s "~{bam}" "~{filename}"
    
    # 2. Index the local symlink
    samtools index "~{filename}"
  >>>

  runtime {
    docker: docker
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible
  }

  output {
    # samtools index creates <filename>.bam.bai
    File bai = "~{filename}.bai"
  }
}

# TASK 2: Split BAM by Read Group
task SplitBamByRG {
  input {
    File bam
    File bai
    String docker
    Int disk_gb
    Int memory_gb
    Int cpu
    Int preemptible
  }

  String filename = basename(bam)

  command <<<
    set -e
    # 1. Symlink BAM and BAI to current directory
    # This forces them to exist side-by-side, which samtools requires
    ln -s "~{bam}" "~{filename}"
    ln -s "~{bai}" "~{filename}.bai"

    # 2. Split using the local symlink
    # -f: Naming format (%* = original basename, %# = Read Group Index)
    # -u: Unaccounted reads (if any) go here
    samtools split -f '%*_%#.bam' -u unassigned.bam "~{filename}"
  >>>

  runtime {
    docker: docker
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible
  }

  output {
    # Capture all BAMs generated by the split
    Array[File] split_bams = glob("*_*.bam")
  }
}

# TASK 3: Convert BAM to FASTQ
task BamToFastq {
  input {
    File bam
    String docker
    Int disk_gb
    Int memory_gb
    Int cpu
    Int preemptible
  }

  # Extract the filename base (e.g., "sample_RG1.bam" -> "sample_RG1")
  String basename = basename(bam, ".bam")

  command <<<
    set -e
    # PIPELINE EXPLANATION:
    # 1. collate: Shuffles reads so pairs are together (required for fastq)
    #    -u: Uncompressed output (faster for pipe)
    #    -O: Output to stdout
    # 2. fastq: Converts to FASTQ format
    #    -1/-2: Paired output files
    #    -0/-s: Discard singletons/orphans (GDC recommendation)
    #    -n:    Use 'standard' /1 /2 suffixes in header
    #    -O:    RESTORE ORIGINAL QUALITIES (Critical GDC Requirement)
    
    samtools collate -u -O ~{bam} | \
    samtools fastq \
      -1 ~{basename}_R1.fastq.gz \
      -2 ~{basename}_R2.fastq.gz \
      -0 /dev/null \
      -s /dev/null \
      -n \
      -O \
      -
  >>>

  runtime {
    docker: docker
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible
  }

  output {
    File r1 = "~{basename}_R1.fastq.gz"
    File r2 = "~{basename}_R2.fastq.gz"
  }
}

# TASK 4: Merge Fastq files
task MergeFastqs {
  input {
    Array[File] r1_files
    Array[File] r2_files
    String base_name
    String docker
    Int disk_gb
    Int memory_gb
    Int cpu
    Int preemptible
  }

  command <<<
    set -e
    # Concatenate the gzipped files (valid for bgzip/gzip)
    cat ~{sep=" " r1_files} > "~{base_name}_R1.fastq.gz"
    cat ~{sep=" " r2_files} > "~{base_name}_R2.fastq.gz"
  >>>

  runtime {
    docker: docker
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible
  }

  output {
    File merged_r1 = "~{base_name}_R1.fastq.gz"
    File merged_r2 = "~{base_name}_R2.fastq.gz"
  }
}

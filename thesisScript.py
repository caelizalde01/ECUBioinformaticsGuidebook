
import argparse
from ast import If
import json
import subprocess
import os
import sys
import glob
import pandas as pd
import zipfile
import shutil
import datetime
import time


if not os.environ.get("_RPY2_ENV_FIXED"):
    r_home = subprocess.run(["R", "RHOME"], capture_output=True, text=True).stdout.strip()
    if r_home:
        r_lib = os.path.join(r_home, "lib")
        ld_dirs = [d for d in os.environ.get("LD_LIBRARY_PATH", "").split(":") if d]
        os.environ["R_HOME"] = r_home
        os.environ["LD_LIBRARY_PATH"] = ":".join([r_lib] + [d for d in ld_dirs if d != r_lib])
        os.environ["_RPY2_ENV_FIXED"] = "1"
        os.execve(sys.executable, [sys.executable] + sys.argv, os.environ)

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr, data
import rpy2.robjects.packages as rpackages



def run_command(command, description):
    """
    Helper function to run shell commands with error handling.
    """
    print(f"\n[STATUS] Starting: {description}")
    print(f"Command: {' '.join(command)}")
    try:
        subprocess.run(command, check=True)
        print(f"[SUCCESS] Completed: {description}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed during: {description}")
        print(f"Error details: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print(f"[ERROR] Executable not found for command: {command[0]}")
        print("Please ensure the tool is installed and in your PATH.")
        sys.exit(1)




def check_tools(tools):
    """
    Verify that necessary tools are installed and available in PATH.
    """
    missing = []
    for tool in tools:
        if shutil.which(tool) is None:
            missing.append(tool)
    
    if missing:
        print(f"[ERROR] The following tools are missing from your PATH: {', '.join(missing)}")
        print("Please install them (e.g., via Conda) before running this script.")
        sys.exit(1)



def find_fastq_pairs(directory):
    """
    Scans a directory for paired fastq files matching the pattern *_1.fastq* and *_2.fastq*.
    Returns a list of tuples: [(sample_name, r1_path, r2_path), ...]
    """
    pairs = []
    # Search for forward reads (supports .fastq, .fq, .fastq.gz, .fq.gz)
    # Using a broad pattern to catch the _1 marker
    forward_files = glob.glob(os.path.join(directory, "*_1.fastq*")) + \
                    glob.glob(os.path.join(directory, "*_1.fq*"))
    
    for r1 in forward_files:
        # Determine the expected reverse file name
        # This replaces the last occurrence of _1 with _2
        base_name = os.path.basename(r1)
        
        # Simple string replacement for the pairing suffix
        if "_1.fastq" in r1:
            r2 = r1.replace("_1.fastq", "_2.fastq")
            sample_id = base_name.split("_1.fastq")[0]
        elif "_1.fq" in r1:
            r2 = r1.replace("_1.fq", "_2.fq")
            sample_id = base_name.split("_1.fq")[0]
        else:
            continue # specific extension not matched logic above, skip

        if os.path.exists(r2):
            pairs.append((sample_id, r1, r2))
        else:
            print(f"[WARNING] Found forward read {base_name} but no matching reverse read found at {r2}. Skipping.")
    
    return sorted(pairs)








def main():

    script_start_time = datetime.datetime.now()


    print(f"Script Start Time: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}")

    parser = argparse.ArgumentParser(description="Automated RNA-Seq Pipeline: Download -> Index -> Align -> Count -> QC -> R Analysis")
                                
    # Required arguments
    parser.add_argument("-a", required=True, help="NCBI RefSeq Assembly Accession (e.g., GCF_011100685.1 for Canis lupus familiaris)")
    parser.add_argument("-fq", required=True, help="Directory containing paired FASTQ files (e.g., sample_1.fastq.gz)")



    # Optional arguments
    parser.add_argument("--threads", default="10", help="Number of threads to use (default: 10)")
    parser.add_argument("--output_dir", default="rnaseq_output", help="Directory for output files")
    parser.add_argument("--keep_zip", action="store_true", help="Keep the downloaded NCBI zip file")
    
    args = parser.parse_args()




    # 0. Pre-flight check
    preflight_start_time = datetime.datetime.now()
    print(f"Preflight Start: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S)")}")
    print(f"Run Time: {preflight_start_time - script_start_time}\n\n")

    required_tools = ["datasets", "hisat2-build", "hisat2", "featureCounts", "multiqc", "fastqc", "rpy2", "r"]
    check_tools(required_tools)

    preflight_end_time = datetime.datetime.now()
    print(f"Preflight End: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S)")}")
    print(f"Run Time: {preflight_end_time - script_start_time}\n\n")


    # 1. Verify FASTQ inputs
    ver_fastq_start_time = datetime.datetime.now()
    print(f"fq Verification Start: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S)")}")
    print(f"Run Time: {ver_fastq_start_time - script_start_time}\n\n")

    if not os.path.isdir(args.fq):
        print(f"[ERROR] FASTQ directory not found: {args.fq}")
        sys.exit(1)
        
    sample_pairs = find_fastq_pairs(args.fq)
    if not sample_pairs:
        print(f"[ERROR] No paired FASTQ files found in {args.fq} matching pattern *_1.fastq*/_2.fastq*")
        sys.exit(1)
    
    print(f"[STATUS] Found {len(sample_pairs)} sample pairs to process:")
    for s, r1, r2 in sample_pairs:
        print(f"  - {s}")

    ver_fastq_end_time = datetime.datetime.now()
    print(f"fq Verification End: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S)")}")
    print(f"Run Time: {ver_fastq_end_time - script_start_time}\n\n")


    # Setup directories
    setup_dir_start_time = datetime.datetime.now()
    print(f"Setup Dir Start: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S)")}")
    print(f"Run Time: {setup_dir_start_time - script_start_time}\n\n")

    os.makedirs(args.output_dir, exist_ok=True)
    genome_dir = os.path.join(args.output_dir, "genome_data")
    counts_dir = os.path.join(args.output_dir, "counts")
    align_dir = os.path.join(args.output_dir, "alignments")
    qc_dir = os.path.join(args.output_dir, "qc")
    multiqc_dir = os.path.join(args.output_dir, "qc", "multiqc_report")
    fastqc_dir = os.path.join(args.output_dir, "qc", "fastqc_report")
    os.makedirs(qc_dir, exist_ok=True)
    os.makedirs(fastqc_dir, exist_ok=True)
    os.makedirs(multiqc_dir, exist_ok=True)
    os.makedirs(counts_dir, exist_ok=True)
    os.makedirs(align_dir, exist_ok=True)

    


    setup_dir_end_time = datetime.datetime.now()
    print(f"Setup Dir End: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S)")}")
    print(f"Run Time: {setup_dir_end_time - script_start_time}\n\n")


    # 2. Download Genome Data
    genome_download_start_time = datetime.datetime.now()
    print(f"Genome Download Start: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S)")}")
    print(f"Run Time: {genome_download_start_time - script_start_time}\n\n")
    zip_filename = os.path.join(args.output_dir, "ncbi_dataset.zip")
    
    # Only download if we don't already have the extracted data
    if not os.path.exists(genome_dir):
        download_cmd = [
            "datasets", "download", "genome", "accession", args.a,
            "--include", "gtf,genome", # genome implies fasta
            "--filename", zip_filename
        ]
        run_command(download_cmd, "Downloading RefSeq Genome and GTF")

        print(f"[STATUS] Unzipping {zip_filename}...")
        try:
            with zipfile.ZipFile(zip_filename, 'r') as zip_ref:
                zip_ref.extractall(genome_dir)
        except zipfile.BadZipFile:
            print("[ERROR] Downloaded file is not a valid zip file.")
            sys.exit(1)
        
        if not args.keep_zip:
            os.remove(zip_filename)
    else:
        print(f"[STATUS] Genome directory {genome_dir} exists. Skipping download.")


    genome_download_end_time = datetime.datetime.now()
    print(f"Genome Download End: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S)")}")
    print(f"Run Time: {genome_download_end_time - script_start_time}\n\n")



    # Locate FASTA and GTF files
    locate_genome_files_start_time = datetime.datetime.now()
    print(f"Locate Genome Files Start: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S)")}")
    print(f"Run Time: {locate_genome_files_start_time - script_start_time}\n\n")

    search_path = os.path.join(genome_dir, "ncbi_dataset", "data", args.a)
    fasta_files = glob.glob(os.path.join(search_path, "*.fna")) + glob.glob(os.path.join(search_path, "*.fasta"))
    gtf_files = glob.glob(os.path.join(search_path, "*.gtf"))

    if not fasta_files:
        print(f"[ERROR] No FASTA (.fna/.fasta) file found in {search_path}")
        sys.exit(1)
    if not gtf_files:
        print(f"[ERROR] No GTF file found in {search_path}")
        sys.exit(1)

    ref_fasta = fasta_files[0]
    ref_gtf = gtf_files[0]
    
    print(f"Found FASTA: {ref_fasta}")
    print(f"Found GTF: {ref_gtf}")

    locate_genome_files_end_time = datetime.datetime.now()
    print(f"Locate Genome Files End: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}")
    print(f"Run Time: {locate_genome_files_end_time - script_start_time}\n\n")

    # 3. Build HISAT2 Index
    hisat2_index_build_start_time = datetime.datetime.now()
    print(f"HISAT2 Index Build Start: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S)")}")
    print(f"Run Time: {hisat2_index_build_start_time - script_start_time}\n\n")

    index_prefix = os.path.join(genome_dir, f"{args.a}_index")
    if not os.path.exists(f"{index_prefix}.1.ht2"):
        build_cmd = [
            "hisat2-build",
            "-p", args.threads,
            ref_fasta,
            index_prefix
        ]
        run_command(build_cmd, "Building HISAT2 Index")
    else:
        print("[STATUS] HISAT2 index found. Skipping build.")
    hisat2_index_build_end_time = datetime.datetime.now()
    print(f"HISAT2 Index Build End: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S)")}")
    print(f"Run Time: {hisat2_index_build_end_time - script_start_time}\n\n")


    # 3A. Run FastQC on all raw reads

    FastQC_Start_time = datetime.datetime.now()

    print("\n" + "="*30)
    print(f"STARTING FASTQC (time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')})")
    print("="*30)
    
    # Extract all R1 and R2 file paths from your sample_pairs list
    all_fastq_files = []
    for s, r1, r2 in sample_pairs:
        all_fastq_files.extend([r1, r2])
        
    fastqc_cmd = [
        "fastqc",
        "-t", args.threads,  # Use the threads argument
        "-o", fastqc_dir     # Output to the folder we just created
    ] + all_fastq_files      # Append all the sequence files to the command
    
    run_command(fastqc_cmd, "Running FastQC on raw reads")

    FastQC_End_time = datetime.datetime.now()


    # 4. Process Loop
    print("\n" + "="*30)
    print(f"STARTING BATCH PROCESSING (time: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")})")
    print("="*30)

    # Prepare lists for featureCounts
    sam_list = []
    counts_list = []
    id_number = 1


    for sample_id, r1, r2 in sample_pairs:
        print(f"\n---> Processing Sample: {sample_id} ({id_number} of {len(sample_pairs)})\nStart Time: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}")

         # Define output file paths
        
        sam_output = os.path.join(align_dir, f"{sample_id}.sam")
        sam_list.append(sam_output)

        counts_output = os.path.join(counts_dir, f"{sample_id}_counts.tsv")
        counts_list.append(counts_output)

        # A. Run HISAT2 Alignment

        summary_log = os.path.join(align_dir, f"{sample_id}_hisat2_summary.txt") # Create log path

        align_cmd = [
            "hisat2",
            "-p", args.threads,
            "-x", index_prefix,
            "-1", r1,
            "-2", r2,
            "-S", sam_output,
            "--summary-file", summary_log
        ]
        run_command(align_cmd, f"Aligning {sample_id}")

        print(f"Completed {sample_id} ({id_number}/{len(sample_pairs)}): {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}")
        id_number += 1


        

    # B. Run featureCounts
    fc_start_time = datetime.datetime.now()

    print(f"\n{('='*30)}\nAbout to attempt featureCounts for all samples again...\nCross thy fingers...\n")
    
    print(f"featureCounts Start Time: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n{('='*30)}")
    print(f"Run Time: {ver_fastq_start_time - script_start_time}\n\n")


    base_count_cmd = [
        "featureCounts",
        "-T", args.threads,
        "-p",             # Paired-end
        #"--countReadPairs", # Count fragments instead of reads --- commented out for testing
        #"-t", "exon",     # Feature type --- Commented out for testing
        #"-g", "gene_id",  # Attribute type --- Commented out for testing
        "-a", ref_gtf,
        "-o", "rnaseq_output/counts/allignedReads.tsv", # Output file
        ]
    
    count_cmd = base_count_cmd + sam_list

    run_command(count_cmd, f"Counting {len(sam_list)} samples with featureCounts")

    fc_end_time = datetime.datetime.now()
    print(f"featureCounts End Time: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S)")}")

    # C. Run MultiQC
    print(f"\n{('='*30)}\nStarting MultiQC Analysis...\n")
    qc_start_time = datetime.datetime.now()
    print(f"MultiQC Start Time: {qc_start_time.strftime("%Y-%m-%d %H:%M:%S")}\n{('='*30)}")
    

    qc_cmd = [
        "multiqc",
        "-o", multiqc_dir,
        "."
    ]

    run_command(qc_cmd, "Running MultiQC on featureCounts output")

    qc_end_time = datetime.datetime.now()
    print(f"MultiQC End Time: {qc_end_time.strftime("%Y-%m-%d %H:%M:%S")}\n{('='*30)}")

    # Do R analysis
    r_start_time = datetime.datetime.now()

    print("\n" + "="*30)
    print("Starting R Analysis")
    print("="*30 + "\n")

    base = importr('base')
    deseq = importr('DESeq2')
    pheatmap = importr('pheatmap')
    tidyverse = importr('tidyverse')
    deseq2 = importr('DESeq2')
    ggplot2 = importr('ggplot2')
    pheatmap = importr('pheatmap')
    dplyr = importr('dplyr')
    fgsea = importr('fgsea')
    msigdbr = importr('msigdbr')
    lazyeval = importr('lazyeval')
    EnhancedVolcano = importr('EnhancedVolcano')

    counts_file = "rnaseq_output/counts/allignedReads.tsv"
    metadata_file = "inputFiles/metadata.csv"

    rDiag_dir = "rnaseq_output/rFigures/diagnostics"
    os.makedirs(rDiag_dir, exist_ok=True)


        # 1. Read the metadata to identify all unique conditions
    metadata = pd.read_csv(metadata_file)


    unique_conditions = metadata['Condition'].unique()

    # 2. Define your control group and extract the experimental groups
    control_group = "Control"
    experimental_groups = [cond for cond in unique_conditions if cond != control_group]



    print("prepping data in R...")

    r_dataPrep = rf"""
    library(DESeq2)

    # 1. Read the featureCounts file
    raw_counts <- read.table("{counts_file}", header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 1)

    # 2. Extract the count matrix (columns 7 to the end)
    count_matrix <- raw_counts[, 7:ncol(raw_counts)]
    rownames(count_matrix) <- raw_counts$Geneid

    # 3. Clean up the sample names in the column headers
    clean_sample_names <- colnames(count_matrix)
    clean_sample_names <- gsub("^.*\\.alignments\\.", "", clean_sample_names) 
    clean_sample_names <- gsub("\\.sam$", "", clean_sample_names)            
    colnames(count_matrix) <- clean_sample_names

    # 4. Load the real metadata table directly 
    metadata <- read.csv("{metadata_file}", stringsAsFactors = FALSE)
    rownames(metadata) <- metadata$SampleID
    metadata <- metadata[colnames(count_matrix), ] # Ensure matching order
    metadata$Condition <- as.factor(metadata$Condition)

    # 5. Create DESeq2 Object
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~ Condition)

    # Run DESeq to calculate dispersions
    dds <- DESeq(dds)

    # 6. Variance Stabilizing Transformation 
    vsd <- vst(dds, blind = FALSE)
    """

    # Execute the code in R
    robjects.r(r_dataPrep)

    print("data prepped...")



    print("graphing pca in R...")

    r_pca = rf"""
    plotPCA(vsd, intgroup = "Condition") + 
    theme_minimal() + 
    ggtitle("PCA of RNA-seq Samples") +
    geom_text(aes(label=name), size=3)

    ggsave("pcaGraph.png", path = "/home/caden/thesisSandbox/script/rnaseq_output/rFigures/diagnostics")
    """

    # Execute the code in R
    robjects.r(r_pca)

    print("pca graphed...")



    r_dispersion = rf"""


    # 1. Open a PNG device and specify the full path and filename
    png("/home/caden/thesisSandbox/script/rnaseq_output/rFigures/diagnostics/dispersionEstimates.png", width=800, height=600)
    
    # 2. Draw the base R plot (it routes directly to the PNG, skipping the pop-up)
    plotDispEsts(dds, main="Dispersion Estimates")
    
    # 3. Close the device to finalize and save the file
    dev.off()
    """

    # Execute the code in R
    robjects.r(r_dispersion)

    print("dispersion graph made...")




    # Beginning of the loop for DEGs, Volcano, and GSEA







    # 1. Automatically extract the species from the NCBI assembly report
    jsonl_path = 'rnaseq_output/genome_data/ncbi_dataset/data/assembly_data_report.jsonl' # Update with your path
    with open(jsonl_path, 'r') as f:
        # A .jsonl file contains one JSON object per line. We just need the first line.
        assembly_data = json.loads(f.readline())

    # Extract the species name (e.g., "Canis lupus familiaris")
    target_species = assembly_data['organism']['organismName']
    print(f"Target species extracted from assembly report: {target_species}")

    # 2. Read the metadata to identify conditions (No need to change your metadata.csv!)


    # Optional: Validate the extracted species against msigdbr to catch spelling mismatches early
    try:
        robjects.r(f"""
        library(msigdbr)
        valid_species <- msigdbr_species()$species_name
        if (!("{target_species}" %in% valid_species)) {{
            stop("Error: '{target_species}' is not a valid species in msigdbr.")
        }}
        """)
        print("Species validation passed!")
    except Exception as e:
        print(f"Species validation failed: {e}")
        # You could uncomment the line below to stop the script if validation fails
        exit(1)



    # Make exp group directories

    for group in experimental_groups:
        gsea_dir = f"/home/caden/thesisSandbox/script/rnaseq_output/rFigures/results/{group}/gsea"
        deg_dir = f"/home/caden/thesisSandbox/script/rnaseq_output/rFigures/results/{group}/deg"
        volcano_dir = f"/home/caden/thesisSandbox/script/rnaseq_output/rFigures/results/{group}/volcanoPlot"
        

        os.makedirs(gsea_dir, exist_ok=True)
        os.makedirs(deg_dir, exist_ok=True)
        os.makedirs(volcano_dir, exist_ok=True)






        


    # 3. Loop through each experimental group and run the DESeq2 extraction
    for exp_group in experimental_groups:

        print(f"Starting analysis on {exp_group} vs {control_group}...")


        
        print(f"Generating DEGs for {exp_group} vs {control_group}...")
        
        # Create the R script string dynamically using an f-string
        # Notice we insert {exp_group} and {control_group} into the contrast and filename
        r_difExp = f"""
        # 1. Extract the results for your specific comparison
        res <- results(dds, contrast=c("Condition", "{exp_group}", "{control_group}"))

        # 2. Order the results by the adjusted p-value
        res <- res[order(res$padj), ]

        # 3. Filter for significance (FDR < 0.05) and effect size (absolute log2FC > 1)
        sig_res <- subset(na.omit(res), padj < 0.05 & abs(log2FoldChange) > 1)

        # 4. Save the significant genes to a CSV file for your records
        write.csv(as.data.frame(sig_res), file="/home/caden/thesisSandbox/script/rnaseq_output/rFigures/results/{exp_group}/deg/Significant_DEGs_{exp_group}_vs_{control_group}.csv")
        """

        robjects.r(r_difExp)
        print(f"Completed DEGs for {exp_group} vs {control_group}!")




        
        print(f"Making Volcano Plot for {exp_group} vs {control_group}...")
        
        r_volcano = f"""
        top5 <- head(rownames(subset(na.omit(as.data.frame(res)), abs(log2FoldChange) > 1 & padj < 0.05)), 5)

        # Generate the plot
        EnhancedVolcano(res,
            lab = rownames(res),
            selectLab = top5,            # <-- only label these genes
            x = 'log2FoldChange',
            y = 'padj',
            title = '{exp_group} vs {control_group}',
            legendLabels = c("Not Influential & Insignificant", "Influential", "Significant", "Influential & Significant"),
            caption = NULL,
            pCutoff = 0.05,
            FCcutoff = 1,
            pointSize = 2.0,
            labSize = 3.0,
            drawConnectors = TRUE,       # draws lines from label to point
            widthConnectors = 0.5,
            arrowheads = FALSE,
            min.segment.length = .1,
            col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
            legendPosition = 'right')

        ggsave("volcanoPlot-{exp_group}.png", path = "/home/caden/thesisSandbox/script/rnaseq_output/rFigures/results/{exp_group}/volcanoPlot")
        """

        robjects.r(r_volcano)
        print(f"Completed volcano plot for {exp_group} vs {control_group}!")


        print(f"making GSEA plot for {exp_group} vs {control_group}...")




        r_gsea = f"""
        library(fgsea)
        library(dplyr)
        library(ggplot2)
        library(msigdbr)

        # ==========================================
        # 1. Prepare the Ranked Gene List
        # ==========================================
        res_clean <- res[!is.na(res$stat), ]
        ranks <- res_clean$stat
        ranks <- ranks + rnorm(length(ranks), sd=1e-10) # Jitter to prevent ties
        names(ranks) <- rownames(res_clean)
        ranks <- sort(ranks, decreasing = TRUE)

        # ==========================================
        # 2. Load Gene Sets (Pathways) using msigdbr
        # ==========================================
        pathways_df <- msigdbr(species = "{target_species}", collection = "H")
        pathways <- split(x = pathways_df$gene_symbol, f = pathways_df$gs_name)

        # ==========================================
        # 3. Run fgsea
        # ==========================================
        fgseaRes <- fgsea(pathways = pathways, 
                        stats = ranks,
                        minSize = 15,
                        maxSize = 500)
        
        fgseaRes <- fgseaRes[order(fgseaRes$NES, decreasing = TRUE), ]
        
        # Flatten the leadingEdge list column to avoid the 'EncodeElement' CSV error
        fgseaRes_save <- fgseaRes
        fgseaRes_save$leadingEdge <- sapply(fgseaRes_save$leadingEdge, paste, collapse = ";")
        
        write.csv(as.data.frame(fgseaRes_save), 
                file="/home/caden/thesisSandbox/script/rnaseq_output/rFigures/results/{exp_group}/gsea/GSEA_{exp_group}_vs_{control_group}2.csv", 
                row.names=FALSE)

                
        
        # ==========================================
        # 4. Visualize and Save the Results
        # ==========================================
                
        if(nrow(fgseaRes[ES > 0]) > 0) {{
            topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=1), pathway]
            topPathwayRes <- fgseaRes[pathway == topPathwaysUp, ]

            p <- plotEnrichment(pathways[[topPathwaysUp]], ranks) + 
                labs(
                    title = paste("Enrichment Plot:", topPathwaysUp),
                    subtitle = paste0(
                        "p-value = ", signif(topPathwayRes$pval, 3),
                        "  |  padj = ", signif(topPathwayRes$padj, 3),
                        "  |  NES = ", round(topPathwayRes$NES, 2)
                        )
                    ) +
                theme_minimal()

                ggsave(filename="GSEA_{exp_group}_vs_{control_group}2.png", plot=p, 
                path="/home/caden/thesisSandbox/script/rnaseq_output/rFigures/results/{exp_group}/gsea", 
                width=8, height=6, bg="white")
        }} else {{
            print("No pathways had an Enrichment Score > 0 to plot.")
        }}

        """


        robjects.r(r_gsea)
        print(f"Completed GSEA plot for {exp_group} vs {control_group}!")

        print(f"{exp_group} analysis complete...")

    r_end_time = datetime.datetime.now()

    print(f"\n{('='*30)}\nR Analysis Complete\n Total time: {r_end_time - r_start_time}\n {('='*30)}\n ")



    print(f"\n{('='*30)}\nJob's done!\n{('='*30)}\n ")
    print("                             ..                     \n                             +#+:                                                -#%#=\n                             .*##*-                 \n                             .+##*+                 \n                             -+***=                 \n                            :++**+-                 \n                           :**+*+-                  \n                          :****+=.                  \n                         =#**+=-.                   \n                       .*#######*****+=-:.          \n                      :#####@@@@%%%%%%%%##*=        \n                     =##**##%%##****#######*:       \n                   :+****##*##%##*++++++++=:        \n                  -********##%@@@%#######*=:        \n                :=+********###%%%#****#****+:       \n              :**+*****##*******+++==+++===-        \n    ..::--===*##****###****++**##**++++==-:         \n**###%%%%%#####*####**********#%@%######**+:        \n###%%#####%######**************###*++=++==-.        \n##########*##***++++*******++++**##*=+++==          \n##########****++++++++*******+**#%%#****+=.         \n#######%%%%###**++**+++++*******##*+++==-.          \n########%%%####**+++++++++++++++++==-::             \n**++++++******#**++++++++++=====--:.                \n+========---===+++=========-::.                     \n=-::..")
    script_end_time = datetime.datetime.now()

    print("\n" + "="*30)
    print("PIPELINE COMPLETE")
    print(f"Counts files location: {counts_dir}")
    print(f"Alignment files location: {align_dir}")
    print("="*30 + "\n")

    print("="*30)
    print(f"Summary of Pipeline\n")
    print(f"Total Time:\n {script_start_time - script_end_time}\n ({script_start_time.strftime("%Y-%m-%d %H:%M:%S")} - {script_end_time.strftime("%Y-%m-%d %H:%M:%S")})\n")
    print(f"Pre-flight Check Time:\n {preflight_start_time - preflight_end_time}\n ({preflight_start_time.strftime("%Y-%m-%d %H:%M:%S")} - {preflight_end_time.strftime("%Y-%m-%d %H:%M:%S")})\n")
    print(f"FASTQ Verification Time:\n {ver_fastq_start_time - ver_fastq_end_time}\n ({ver_fastq_start_time.strftime("%Y-%m-%d %H:%M:%S")} - {ver_fastq_end_time.strftime("%Y-%m-%d %H:%M:%S")})\n")
    print(f"Directory Setup Time:\n {setup_dir_start_time - setup_dir_end_time}\n ({setup_dir_start_time.strftime("%Y-%m-%d %H:%M:%S")} - {setup_dir_end_time.strftime("%Y-%m-%d %H:%M:%S")})\n")
    print(f"Genome Download Time:\n {genome_download_start_time - genome_download_end_time}\n ({genome_download_start_time.strftime("%Y-%m-%d %H:%M:%S")} - {genome_download_end_time.strftime("%Y-%m-%d %H:%M:%S")})\n")
    print(f"Locate Genome Files Time:\n {locate_genome_files_start_time - locate_genome_files_end_time}\n ({locate_genome_files_start_time.strftime("%Y-%m-%d %H:%M:%S")} - {locate_genome_files_end_time.strftime("%Y-%m-%d %H:%M:%S")})\n")
    
    print(f"FastQC Start Time:\n {FastQC_Start_time - FastQC_End_time}\n ({FastQC_Start_time.strftime("%Y-%m-%d %H:%M:%S")} - {FastQC_End_time.strftime("%Y-%m-%d %H:%M:%S")})\n")

    print(f"HISAT2 Index Build Time:\n {hisat2_index_build_start_time - hisat2_index_build_end_time}\n ({hisat2_index_build_start_time.strftime("%Y-%m-%d %H:%M:%S")} - {hisat2_index_build_end_time.strftime("%Y-%m-%d %H:%M:%S")})\n")
    print(f"featureCounts Time:\n {fc_start_time - fc_end_time}\n ({fc_start_time.strftime("%Y-%m-%d %H:%M:%S")} - {fc_end_time.strftime("%Y-%m-%d %H:%M:%S")})\n")
    print(f"MultiQC Time:\n {qc_start_time - qc_end_time}\n ({qc_start_time.strftime("%Y-%m-%d %H:%M:%S")} - {qc_end_time.strftime("%Y-%m-%d %H:%M:%S")})\n")
    print(f"R Analysis Time:\n {r_start_time - r_end_time}\n ({r_start_time.strftime("%Y-%m-%d %H:%M:%S")} - {r_end_time.strftime("%Y-%m-%d %H:%M:%S")})\n")
    print("="*30)



if __name__ == "__main__":
    main()


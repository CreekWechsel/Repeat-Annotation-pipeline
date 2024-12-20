TE_Annotation_Pipeline
Introduction
TE_Annotation_Pipeline.py is an integrated pipeline that annotates Transposable Elements (TEs) in a given genome sequence. It incorporates multiple steps and external tools for LTR prediction, repeat de novo prediction, homology-based filtering, repeat library merging, final annotation, masking of TE regions, and simple repeat detection. This pipeline provides a one-stop solution for repetitive sequence annotation workflows.

Features and Steps
The pipeline consists of the following nine steps:

LTRharvest for LTR Prediction (01.ltrharvest):
Use gt ltrharvest to perform an initial LTR-like transposable element prediction on the genome.

LTRfinder for LTR Prediction (02.ltrfinder):
Use LTR_FINDER_parallel to predict LTRs, complementing the results from LTRharvest.

LTR_retriever for Filtering and Integration (03.LTR_retriever):
Integrate and filter results from the first two steps using LTR_retriever to obtain a high-quality LTR library.

RepeatModeler de novo Prediction (04.RepeatModeler):
Use RepeatModeler on the partially masked genome to perform de novo prediction of repetitive elements, generating a repeat library (Modeler.lib).

Extracting Homologous Repeats with famdb.py (Within 05.RepeatMasker, Step 5):
Use famdb.py to extract related repeat sequences based on species information.

Merging Three LIBs (Within 05.RepeatMasker, Step 6):
Merge the RepeatModeler library, the famdb.py-extracted library, and the LTR library obtained earlier into a comprehensive repeat library (allRepeats.lib).

Final Repeat Annotation with RepeatMasker (Within 05.RepeatMasker, Step 7):
Use the merged allRepeats.lib with RepeatMasker to annotate repeats in the genome.

Hard and Soft Masking of TE Regions (Within 05.RepeatMasker, Step 8):
Use bedtools maskfasta to apply both hard and soft masking to the TE regions, producing genome.hardmasked.fa and genome.softmasked.fa.

TRF for Simple Repeat Detection (06.trf, Step 9):
Use trf on the genome.hardmasked.fa generated from Step 8 to detect simple tandem repeats, producing .dat and .gff files.

After each step completes, a .stepN_done checkpoint file is created in the corresponding directory. Using the --resume option, if these checkpoint files are found, the pipeline will skip the completed steps and resume from the last unfinished step.

External Tool Dependencies
This pipeline depends on several external tools that must be installed and available in your PATH:

gt (genometools)
LTR_FINDER_parallel
LTR_retriever
RepeatModeler
RepeatMasker
famdb.py (included with EDTA or other packages)
bedtools
trf
The script checks for these tools at runtime and will exit with an error if any are missing.

Parameters
bash
复制代码
python3 TE_Annotation_Pipeline.py --genome <genome.fasta> --sp <species_name> [options]
Required Arguments:

--genome <str>: Path to the input genome FASTA file
--sp <str>: Species name used for famdb.py filtering
Optional Arguments (with defaults):

--threads <int>: Number of threads (default: 32)
--minlenltr <int>: Minimum LTR length for LTRharvest (default: 100)
--maxlenltr <int>: Maximum LTR length for LTRharvest (default: 6000)
--mindistltr <int>: Minimum distance between LTRs (default: 1500)
--maxdistltr <int>: Maximum distance between LTRs (default: 25000)
--similar <int>: Similarity parameter for LTRharvest (default: 85)
--mintsd <int>: Minimum TSD length (default: 5)
--maxtsd <int>: Maximum TSD length (default: 5)
--vic <int>: VIC parameter for LTRharvest (default: 10)
--database <str>: Database name for BuildDatabase and RepeatModeler (default: {sp}_DB)
--resume: Attempt to resume from previously completed steps
Example Usage
Basic run with default parameters and 32 threads:

bash
复制代码
python3 TE_Annotation_Pipeline.py --genome genome.fasta --sp MySpecies
Specify more threads, adjust LTR parameters, and resume after interruption:

bash
复制代码
python3 TE_Annotation_Pipeline.py --genome genome.fasta --sp MySpecies --threads 64 --minlenltr 200 --maxlenltr 5000 --resume
Provide a custom database name:

bash
复制代码
python3 TE_Annotation_Pipeline.py --genome genome.fasta --sp MySpecies --database MyDB
Logs and Intermediate Files
The script logs its operations to TE_Annotation_Pipeline.log and creates subdirectories for each step (e.g., 01.ltrharvest, 02.ltrfinder, 03.LTR_retriever, 04.RepeatModeler, 05.RepeatMasker, 06.trf). Each step's intermediate files and results are stored in their respective directories. Checkpoint files (.stepN_done) indicate completed steps and facilitate resuming the pipeline.

Error Handling and Exit Codes
If any external tool returns a non-zero exit code that is not allowed (except for steps where you explicitly relax this condition, e.g., for trf), the script will raise an exception and terminate. Using --resume, you can restart the pipeline, and it will continue from the last successful checkpoint.

License
You may distribute this code under an open-source license such as MIT, GPL, or BSD. Adapt as needed according to your project's requirements.

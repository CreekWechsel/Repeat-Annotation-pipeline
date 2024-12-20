# -*- coding: utf-8 -*-
# Chenzj
# a604249194@126.com
# SNNU

import os
import sys
import argparse
import logging
import subprocess
import shutil
import glob

def run_cmd(cmd, logger, work_dir=None, allowed_returncodes=None):
    if allowed_returncodes is None:
        allowed_returncodes = [0]
    
    logger.info("Running command: %s", cmd)
    completed = subprocess.run(cmd, shell=True, cwd=work_dir, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if completed.returncode not in allowed_returncodes:
        logger.error("Command failed with return code %s", completed.returncode)
        logger.error("STDOUT: %s", completed.stdout)
        logger.error("STDERR: %s", completed.stderr)
        raise RuntimeError(f"Command failed: {cmd}")
    else:
        logger.debug("Command STDOUT: %s", completed.stdout.strip())
        logger.debug("Command STDERR: %s", completed.stderr.strip())
    return completed.stdout

def ensure_dir(directory, logger):
    if not os.path.exists(directory):
        logger.info("Creating directory: %s", directory)
        os.makedirs(directory, exist_ok=True)

def check_tool_availability(tools, logger):
    logger.info("Checking external tools availability...")
    for tool in tools:
        if shutil.which(tool) is None:
            logger.error("Required tool '%s' not found in PATH.", tool)
            raise FileNotFoundError(f"Tool {tool} not found in PATH")
    logger.info("All required tools are available.")

def create_checkpoint(step_dir, step_num):
    done_file = os.path.join(step_dir, f".step{step_num}_done")
    with open(done_file, 'w') as f:
        f.write("done\n")

def check_checkpoint(step_dir, step_num, resume):
    done_file = os.path.join(step_dir, f".step{step_num}_done")
    if resume and os.path.isfile(done_file):
        return True
    return False

def convert_repeat_to_gff(input_file, output_file, prefix=None, logger=None):
    if logger is None:
        logger = logging.getLogger(__name__)

    def create_id(tag, index):
        return f"{prefix + '_' if prefix else ''}{tag}{index}"

    file_suffix = None
    if input_file.endswith(".dat"):
        file_suffix = "dat"
    elif input_file.endswith(".out"):
        file_suffix = "out"
    elif input_file.endswith(".annot"):
        file_suffix = "annot"
    else:
        raise ValueError("Input file must end with .dat, .out, or .annot")

    line_count = 0
    with open(input_file, 'r') as f:
        for _ in f:
            line_count += 1
    mark_base = line_count + 1

    def parse_dat(file_path):
        chr_name = None
        entries = []
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("Sequence:"):
                    chr_name = line.split()[1]
                    continue
                parts = line.split()
                if len(parts) == 15 and parts[0].isdigit():
                    start, end = parts[0], parts[1]
                    score = parts[7]
                    consensus_size = parts[4]
                    copy_number = parts[3]
                    percent_matches = parts[5]
                    percent_indels = parts[6]
                    consensus_seq = parts[13]
                    repeat_seq = parts[14]
                    entries.append({
                        "chr": chr_name,
                        "source": "TRF",
                        "type": "TandemRepeat",
                        "start": start,
                        "end": end,
                        "score": score,
                        "strand": "+",
                        "attributes": f"ConsensusSize={consensus_size};CopyNumber={copy_number};PercentMatches={percent_matches};PercentIndels={percent_indels};Consensus={consensus_seq};RepeatSeq={repeat_seq}"
                    })
        return entries

    def parse_out(file_path):
        entries = []
        with open(file_path, 'r') as f:
            for line in f:
                line=line.strip()
                parts = line.split()
                if len(parts) < 14:
                    continue
                if not parts[0].isdigit():
                    continue
                te_class = parts[10]
                if ("Low" in te_class) or ("Simple" in te_class) or ("Satellite" in te_class):
                    continue
                score = parts[0]
                chr_name = parts[4]
                start = parts[5]
                end = parts[6]
                strand = "+" if parts[8] == '+' else '-'
                target = parts[9]

                coords = []
                for p in parts[11:]:
                    if '(' not in p and ')' not in p and p.isdigit():
                        coords.append(int(p))
                if len(coords) >= 2:
                    coords.sort()
                    target_start, target_end = coords[0], coords[1]
                else:
                    target_start, target_end = (0,0)

                perc_div = parts[1]
                perc_del = parts[2]
                perc_ins = parts[3]

                attr = f"Target={target} {target_start} {target_end};Class={te_class};PercDiv={perc_div};PercDel={perc_del};PercIns={perc_ins}"
                entries.append({
                    "chr": chr_name,
                    "source": "RepeatMasker",
                    "type": "Transposon",
                    "start": start,
                    "end": end,
                    "score": score,
                    "strand": strand,
                    "attributes": attr
                })
        return entries

    def parse_annot(file_path):
        entries = []
        with open(file_path, 'r') as f:
            for line in f:
                line=line.strip()
                if line.startswith("pValue"):
                    continue
                parts = line.split()
                if len(parts) < 11:
                    continue
                pval = parts[0]
                score = parts[1]
                chr_name = parts[3]
                start = parts[4]
                end = parts[5]
                strand = parts[6]
                target = parts[7]
                te_class = parts[8]
                coords = [int(parts[9]), int(parts[10])]
                coords.sort()
                target_start, target_end = coords[0], coords[1]

                attr = f"Target={target} {target_start} {target_end};Class={te_class};pValue={pval}"
                entries.append({
                    "chr": chr_name,
                    "source": "RepeatProteinMask",
                    "type": "TEprotein",
                    "start": start,
                    "end": end,
                    "score": score,
                    "strand": strand,
                    "attributes": attr
                })
        return entries

    if file_suffix == "dat":
        entries = parse_dat(input_file)
        id_tag = "TR"
    elif file_suffix == "out":
        entries = parse_out(input_file)
        id_tag = "TE"
    elif file_suffix == "annot":
        entries = parse_annot(input_file)
        id_tag = "TP"
    else:
        raise ValueError("Unknown file suffix")

    with open(output_file, 'w') as out_handle:
        mark = mark_base
        for e in entries:
            gff_line = (
                f"{e['chr']}\t{e['source']}\t{e['type']}\t{e['start']}\t{e['end']}"
                f"\t{e['score']}\t{e['strand']}\t.\tID={create_id(id_tag, mark)};{e['attributes']}\n"
            )
            out_handle.write(gff_line)
            mark += 1

    logger.info(f"Converted {input_file} to GFF and wrote to {output_file}.")

def step1_ltrharvest(cwd, genome, indexbase, minlenltr, maxlenltr, mindistltr, maxdistltr, similar,
                     mintsd, maxtsd, vic, logger, resume=False):
    step_dir = "01.ltrharvest"
    step_dir = os.path.join(cwd, step_dir)
    ensure_dir(step_dir, logger)
    if check_checkpoint(step_dir, 1, resume):
        logger.info("Skipping step1 (LTRharvest) as checkpoint found.")
        return

    index_path = os.path.join(step_dir, indexbase)
    out_gff = f"{indexbase}.gff{similar}"
    out_out = f"{indexbase}.out{similar}"
    out_inner = f"{indexbase}.outinner{similar}"
    out_scn = f"{indexbase}.harvest.scn"

    run_cmd(f"gt suffixerator -db {genome} -indexname {index_path} -tis -suf -lcp -des -ssp -sds -dna", logger, work_dir=step_dir)
    run_cmd(
        f"gt ltrharvest -index {index_path} "
        f"-out {out_out} -outinner {out_inner} "
        f"-gff3 {out_gff} "
        f"-minlenltr {minlenltr} -maxlenltr {maxlenltr} "
        f"-mindistltr {mindistltr} -maxdistltr {maxdistltr} "
        f"-mintsd {mintsd} -maxtsd {maxtsd} -vic {vic} "
        f"-similar {similar} > {out_scn}", logger, work_dir=step_dir
    )
    create_checkpoint(step_dir, 1)

def step2_ltrfinder(cwd, genome, threads, logger, resume=False):
    step_dir = os.path.join(cwd, "02.ltrfinder")
    ensure_dir(step_dir, logger)
    if check_checkpoint(step_dir, 2, resume):
        logger.info("Skipping step2 (LTRfinder) as checkpoint found.")
        return

    indexbase = os.path.basename(genome)
    out_scn = f"{indexbase}.finder.combine.scn"
    run_cmd(f"LTR_FINDER_parallel -seq {genome} -threads {threads} -harvest_out -size 1000000 -time 300 -out {out_scn}", logger, work_dir=step_dir)
    create_checkpoint(step_dir, 2)

def step3_ltr_retriever(cwd, genome, logger, threads=32, resume=False):
    step_dir = os.path.join(cwd, "03.LTR_retriever")
    ensure_dir(step_dir, logger)
    if check_checkpoint(step_dir, 3, resume):
        logger.info("Skipping step3 (LTR_retriever) as checkpoint found.")
        return

    indexbase = os.path.basename(genome)
    harvest_scn = os.path.join(cwd, "01.ltrharvest", f"{indexbase}.harvest.scn")
    finder_scn = os.path.join(cwd, "02.ltrfinder", f"{indexbase}.finder.combine.scn")

    run_cmd(f"cat {harvest_scn} {finder_scn} > raw_LTR.scn", logger, work_dir=step_dir)
    run_cmd(f"LTR_retriever -genome {genome} -inharvest raw_LTR.scn -threads {threads} -u 1.01e-8 -out {indexbase}", logger, work_dir=step_dir)
    create_checkpoint(step_dir, 3)

def step4_repeatmodeler(cwd, genome, logger, threads=32, database_name="default_db", resume=False):
    step_dir = os.path.join(cwd, "04.RepeatModeler")
    ensure_dir(step_dir, logger)
    if check_checkpoint(step_dir, 4, resume):
        logger.info("Skipping step4 (RepeatModeler) as checkpoint found.")
        return

    indexbase = os.path.basename(genome)
    ltr_lib_src = os.path.join(cwd, "03.LTR_retriever", f"{indexbase}.LTRlib.fa")
    if not os.path.isfile(ltr_lib_src):
        raise FileNotFoundError(f"Missing LTR library file: {ltr_lib_src}")

    run_cmd(f"cp {ltr_lib_src} LTR.lib", logger, work_dir=step_dir)
    run_cmd(f"RepeatMasker -pa {threads} -no_is -e rmblast -lib LTR.lib -dir 01.RM-1 {genome}", logger, work_dir=step_dir)

    masked_genome = os.path.join(step_dir, "01.RM-1", f"{indexbase}.masked")
    if not os.path.isfile(masked_genome):
        raise FileNotFoundError(f"Expected masked genome file not found: {masked_genome}")

    run_cmd(f"BuildDatabase -name {database_name} -engine rmblast {masked_genome}", logger, work_dir=step_dir)
    run_cmd(f"RepeatModeler -threads {threads} -database {database_name}", logger, work_dir=step_dir)

    rm_dirs = glob.glob(os.path.join(step_dir, "RM_*"))
    consensi_classified = None
    for d in rm_dirs:
        if os.path.isdir(d):
            candidate = os.path.join(d, "consensi.fa.classified")
            if os.path.isfile(candidate):
                consensi_classified = candidate
                break

    if not consensi_classified:
        direct_candidate = os.path.join(step_dir, "consensi.fa.classified")
        if os.path.isfile(direct_candidate):
            consensi_classified = direct_candidate

    if not consensi_classified:
        raise FileNotFoundError("Could not find consensi.fa.classified from RepeatModeler output.")

    run_cmd(f"cp {consensi_classified} Modeler.lib", logger, work_dir=step_dir)
    create_checkpoint(step_dir, 4)

def step5_famdb(cwd, genome, sp, logger, resume=False):
    step_dir = os.path.join(cwd, "05.RepeatMasker")
    ensure_dir(step_dir, logger)
    if check_checkpoint(step_dir, 5, resume):
        logger.info("Skipping step5 (famdb) as checkpoint found.")
        return

    sp_db_path = f"{sp}.db"
    run_cmd(f"famdb.py families --format fasta_name --ancestors --descendants --curated --include-class-in-name {sp} > {sp_db_path}", logger, work_dir=step_dir)
    create_checkpoint(step_dir, 5)

def step6_merge_libs(cwd, genome, sp, logger, resume=False):
    step_dir = os.path.join(cwd, "05.RepeatMasker")
    if check_checkpoint(step_dir, 6, resume):
        logger.info("Skipping step6 (merge libs) as checkpoint found.")
        return

    modeler_lib = os.path.join(cwd, "04.RepeatModeler", "Modeler.lib")
    sp_db = os.path.join(step_dir, f"{sp}.db")
    ltr_lib_path = os.path.join(cwd, "04.RepeatModeler", "LTR.lib")

    for f in [modeler_lib, sp_db, ltr_lib_path]:
        if not os.path.isfile(f):
            raise FileNotFoundError(f"Missing required file: {f}")

    run_cmd(f"cat {modeler_lib} {sp_db} {ltr_lib_path} > allRepeats.lib", logger, work_dir=step_dir)
    create_checkpoint(step_dir, 6)

def step7_repeatmasker_final(cwd, genome, logger, threads=32, resume=False):
    step_dir = os.path.join(cwd, "05.RepeatMasker")
    if check_checkpoint(step_dir, 7, resume):
        logger.info("Skipping step7 (RepeatMasker final) as checkpoint found.")
        return

    local_genome = os.path.join(step_dir, "genome.fasta")
    if not os.path.isfile(local_genome):
        run_cmd(f"cp {genome} genome.fasta", logger, work_dir=step_dir)

    run_cmd(f"RepeatMasker -e rmblast -gff -xsmall -html -norna -source -no_is -pa {threads} -s -lib allRepeats.lib genome.fasta", logger, work_dir=step_dir)

    rm_out = os.path.join(step_dir, "genome.fasta.out")
    if not os.path.isfile(rm_out):
        alt_out = os.path.join(step_dir, "genome.out")
        if os.path.isfile(alt_out):
            rm_out = alt_out
        else:
            raise FileNotFoundError("Cannot find RepeatMasker output file (.out) in 05.RepeatMasker")

    te_repeat_gff = os.path.join(step_dir, "genome.TE_repeat.gff")
    convert_repeat_to_gff(rm_out, te_repeat_gff, prefix=None, logger=logger)

    if not os.path.isfile(te_repeat_gff):
        raise FileNotFoundError("TE repeat GFF file not found after conversion.")

    create_checkpoint(step_dir, 7)

def step8_mask_te(cwd, logger, resume=False):
    step_dir = os.path.join(cwd, "05.RepeatMasker")
    if check_checkpoint(step_dir, 8, resume):
        logger.info("Skipping step8 (TE masking) as checkpoint found.")
        return

    run_cmd("bedtools maskfasta -fi genome.fasta -bed genome.TE_repeat.gff -fo genome.hardmasked.fa", logger, work_dir=step_dir)
    run_cmd("bedtools maskfasta -soft -fi genome.fasta -bed genome.TE_repeat.gff -fo genome.softmasked.fa", logger, work_dir=step_dir)

    create_checkpoint(step_dir, 8)

def step9_trf(cwd, logger, resume=False):
    step_dir = os.path.join(cwd, "06.trf")
    ensure_dir(step_dir, logger)
    if check_checkpoint(step_dir, 9, resume):
        logger.info("Skipping step9 (trf) as checkpoint found.")
        return

    masked_genome = os.path.join(cwd, "05.RepeatMasker", "genome.hardmasked.fa")
    # Allow return codes 0 or 1 to be considered successful
    run_cmd(f"trf {masked_genome} 2 7 7 80 10 50 2000 -h -d", logger, work_dir=step_dir, allowed_returncodes=[0,1])

    trf_dat_1 = os.path.join(step_dir, "genome.hardmasked.fa.2.7.7.80.10.50.2000.dat")
    trf_dat_2 = os.path.join(step_dir, "genome.2.7.7.80.10.50.2000.dat")
    if os.path.isfile(trf_dat_1):
        trf_dat = trf_dat_1
    elif os.path.isfile(trf_dat_2):
        trf_dat = trf_dat_2
    else:
        raise FileNotFoundError("Cannot find TRF output dat file.")

    trf_gff = os.path.join(step_dir, "genome.trf.gff")
    convert_repeat_to_gff(trf_dat, trf_gff, prefix=None, logger=logger)
    create_checkpoint(step_dir, 9)

def main():
    parser = argparse.ArgumentParser(description="A pipeline for TE annotation.")
    parser.add_argument("--genome", required=True, help="Path to the input genome fasta file, e.g., genome.fasta")
    parser.add_argument("--sp", required=True, help="Species name used for famdb.py filtering")
    parser.add_argument("--threads", type=int, default=32, help="Number of threads for parallel tools")
    parser.add_argument("--minlenltr", type=int, default=100, help="Minimum LTR length for LTR_harvest")
    parser.add_argument("--maxlenltr", type=int, default=6000, help="Maximum LTR length for LTR_harvest")
    parser.add_argument("--mindistltr", type=int, default=1500, help="Minimum distance between LTRs")
    parser.add_argument("--maxdistltr", type=int, default=25000, help="Maximum distance between LTRs")
    parser.add_argument("--similar", type=int, default=85, help="Similarity parameter for LTR_harvest")
    parser.add_argument("--mintsd", type=int, default=5, help="Minimum TSD length")
    parser.add_argument("--maxtsd", type=int, default=5, help="Maximum TSD length")
    parser.add_argument("--vic", type=int, default=10, help="VIC parameter for LTR_harvest")
    parser.add_argument("--database", help="Database name for BuildDatabase and RepeatModeler", default=None)
    parser.add_argument("--resume", action="store_true", help="Attempt to resume from previous partial results")

    args = parser.parse_args()

    if args.database is None:
        args.database = f"{args.sp}_DB"

    logger = logging.getLogger("TE_Annotation_Pipeline")
    logger.setLevel(logging.DEBUG)
    
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    fh = logging.FileHandler("TE_Annotation_Pipeline.log", mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    genome = os.path.abspath(args.genome)
    sp = args.sp
    cwd = os.path.abspath(os.getcwd())

    if not os.path.isfile(genome):
        logger.error("Genome file does not exist: %s", genome)
        sys.exit(1)

    required_tools = [
        "gt", "LTR_FINDER_parallel", "LTR_retriever", "RepeatModeler", 
        "RepeatMasker", "famdb.py", "bedtools", "trf"
    ]
    try:
        check_tool_availability(required_tools, logger)
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)

    try:
        step1_ltrharvest(
            cwd, genome, os.path.basename(genome),
            minlenltr=args.minlenltr,
            maxlenltr=args.maxlenltr,
            mindistltr=args.mindistltr,
            maxdistltr=args.maxdistltr,
            similar=args.similar,
            mintsd=args.mintsd,
            maxtsd=args.maxtsd,
            vic=args.vic,
            logger=logger,
            resume=args.resume
        )

        step2_ltrfinder(cwd, genome, threads=args.threads, logger=logger, resume=args.resume)
        step3_ltr_retriever(cwd, genome, threads=args.threads, logger=logger, resume=args.resume)
        step4_repeatmodeler(cwd, genome, threads=args.threads, logger=logger, database_name=args.database, resume=args.resume)
        step5_famdb(cwd, genome, sp, logger=logger, resume=args.resume)
        step6_merge_libs(cwd, genome, sp, logger=logger, resume=args.resume)
        step7_repeatmasker_final(cwd, genome, logger=logger, threads=args.threads, resume=args.resume)
        step8_mask_te(cwd, logger=logger, resume=args.resume)
        step9_trf(cwd, logger=logger, resume=args.resume)

        logger.info("All steps completed successfully.")
    except Exception as e:
        logger.exception("The pipeline encountered an error.")
        sys.exit(1)

if __name__ == "__main__":
    main()

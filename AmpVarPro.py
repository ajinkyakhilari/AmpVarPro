#!/usr/bin/env python3
"""
ONT Pf panel pipeline runner (barcode-aware), no Snakemake required.

Steps per barcode:
  1) Merge FASTQs
  2) fastp filter
  3) minimap2 map -> samtools sort/index
  4) (optional) samtools -L region filter
  5) (optional) iVar trim (amplicon primer trimming)
  6) sort/index (post-trim or region BAM)
  7) per-base coverage cap (<= N) using built-in pysam+numpy logic
  8) samtools coverage (QC)
  9) Clair3
 10) snpEff annotation

Requirements in PATH: fastp, minimap2, samtools, snpEff, run_clair3.sh
Python libs: pysam, numpy

Example:
  python AmpVarPro.py \
    --input-dir /path/fastq_pass \
    --output-dir /path/output \
    --reference-fasta /path/PlasmoDB-61_Pfalciparum3D7_Genome.fasta \
    --regions-file /path/Pf-regions-extract.txt \
    --do-ivar-trim \
    --primer-bed /path/Primer-15.bed \
    --threads 20 \
    --cap-depth 200 \
    --clair3-model-dir /path/models/r941_prom_sup_g5014/ \
    --clair3-preset ont \
    --snpeff-db pf3d7-61
"""

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import List, Dict, Tuple, Optional

# ---------- utilities ----------

def run(cmd: str, log_file: Optional[Path] = None, check: bool = True) -> int:
    """Run a shell command with optional tee to a log file."""
    print(f"[CMD] {cmd}", flush=True)
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with open(log_file, "a") as lf:
            p = subprocess.Popen(cmd, shell=True, stdout=lf, stderr=subprocess.STDOUT, executable="/bin/bash")
            ret = p.wait()
    else:
        ret = subprocess.call(cmd, shell=True, executable="/bin/bash")
    if check and ret != 0:
        raise RuntimeError(f"Command failed ({ret}): {cmd}")
    return ret

def ensure_fai(reference_fasta: Path):
    fai = reference_fasta.with_suffix(reference_fasta.suffix + ".fai")
    if not fai.exists():
        run(f"samtools faidx {shlex.quote(str(reference_fasta))}")

def discover_barcodes(input_dir: Path, pattern: str = "barcode") -> List[Path]:
    bcs = []
    for d in sorted(input_dir.iterdir()):
        if d.is_dir() and d.name.startswith(pattern):
            if any(d.glob("*.fastq.gz")) or any(d.glob("*.fq.gz")):
                bcs.append(d)
    return bcs

def cat_fastqs(srcs: List[Path], dest: Path):
    if not srcs:
        raise ValueError("No FASTQs to merge.")
    files = " ".join(shlex.quote(str(p)) for p in srcs)
    run(f"cat {files} > {shlex.quote(str(dest))}")

# ---------- per-base cap (pysam + numpy) ----------

def cap_bam_per_base(
    in_bam: Path,
    out_bam: Path,
    cap: int,
    regions_bed: Optional[Path] = None,
    seed: Optional[int] = 42,
):
    """
    Cap per-base coverage by skipping whole reads that would push any covered
    position above the cap. Uses uint16 arrays; suitable for caps like 200.
    """
    import numpy as np
    import pysam
    import random
    from collections import defaultdict

    def load_regions_bed(bed_path: Path) -> Dict[str, List[Tuple[int, int]]]:
        regions = defaultdict(list)
        with open(bed_path) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                start = int(parts[1])  # 0-based
                end = int(parts[2])    # end-exclusive
                regions[chrom].append((start, end))
        # merge intervals
        for chrom in list(regions.keys()):
            intervals = sorted(regions[chrom])
            merged = []
            for s, e in intervals:
                if not merged or s > merged[-1][1]:
                    merged.append([s, e])
                else:
                    merged[-1][1] = max(merged[-1][1], e)
            regions[chrom] = [(s, e) for s, e in merged]
        return regions

    def in_any_region(pos: int, intervals: List[Tuple[int, int]]) -> bool:
        lo, hi = 0, len(intervals) - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            s, e = intervals[mid]
            if pos < s:
                hi = mid - 1
            elif pos >= e:
                lo = mid + 1
            else:
                return True
        return False

    def iter_reads_grouped_by_leftpos(bam):
        group = []
        current = None
        for r in bam.fetch(until_eof=False):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            key = (r.reference_id, r.reference_start)
            if current is None:
                current = key
            if key != current:
                yield group
                group = []
                current = key
            group.append(r)
        if group:
            yield group

    if cap <= 0:
        # Just copy with index
        run(f"cp {shlex.quote(str(in_bam))} {shlex.quote(str(out_bam))}")
        run(f"cp {shlex.quote(str(in_bam))}.bai {shlex.quote(str(out_bam))}.bai", check=False)
        return

    if seed is not None:
        import random
        random.seed(seed)

    import pysam
    in_af = pysam.AlignmentFile(str(in_bam), "rb")
    out_af = pysam.AlignmentFile(str(out_bam), "wb", header=in_af.header)

    import numpy as np
    contig_cov = [np.zeros(in_af.get_reference_length(ctg), dtype=np.uint16)
                  for ctg in in_af.references]

    region_map = load_regions_bed(regions_bed) if regions_bed else None

    kept = 0
    skipped = 0

    for batch in iter_reads_grouped_by_leftpos(in_af):
        if seed is not None and len(batch) > 1:
            import random
            random.shuffle(batch)
        for r in batch:
            ref_id = r.reference_id
            ctg = in_af.get_reference_name(ref_id)
            ref_positions = r.get_reference_positions(full_length=False)
            if not ref_positions:
                continue
            if region_map is not None and ctg in region_map:
                intervals = region_map[ctg]
                ref_positions = [p for p in ref_positions if in_any_region(p, intervals)]
                if not ref_positions:
                    out_af.write(r)
                    kept += 1
                    continue
            cov = contig_cov[ref_id]
            if any(cov[p] >= cap for p in ref_positions):
                skipped += 1
                continue
            cov.put(ref_positions, cov.take(ref_positions) + 1)
            out_af.write(r)
            kept += 1

    in_af.close()
    out_af.close()
    pysam.index(str(out_bam))
    print(f"[cap] kept={kept} skipped={skipped} cap={cap}")

# ---------- main per-barcode workflow ----------

def process_barcode(
    bc_dir: Path,
    out_root: Path,
    threads: int,
    ref_fa: Path,
    regions_file: Optional[Path],
    do_ivar: bool,
    primer_bed: Optional[Path],
    fastp_q: int,
    fastp_minlen: int,
    fastp_maxlen: int,
    cap_depth: int,
    cap_seed: int,
    clair3_model_dir: Path,
    clair3_preset: str,
    snpeff_db: str,
):
    bc = bc_dir.name
    bc_out = out_root / bc
    log_dir = bc_out / "logs"
    bc_out.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    # 0) list fastqs
    fastqs = sorted(list(bc_dir.glob("*.fastq.gz")) + list(bc_dir.glob("*.fq.gz")))
    if not fastqs:
        print(f"[WARN] No FASTQs in {bc_dir}, skipping.")
        return

    # 1) merge
    merged_fq = bc_out / f"{bc}.fastq.gz"
    if not merged_fq.exists():
        cat_fastqs(fastqs, merged_fq)

    # 2) fastp
    filt_fq = bc_out / f"{bc}.filt.fastq.gz"
    if not filt_fq.exists():
        run(
            " ".join([
                "fastp",
                f"-i {shlex.quote(str(merged_fq))}",
                f"-o {shlex.quote(str(filt_fq))}",
                f"-q {fastp_q}",
                f"-l {fastp_minlen}",
                f"--length_limit {fastp_maxlen}",
                "-A",
                f"-w {threads}",
                f"--html {shlex.quote(str(bc_out / f'{bc}.fastp.html'))}",
                f"--json {shlex.quote(str(bc_out / f'{bc}.fastp.json'))}",
            ]),
            log_file=log_dir / f"{bc}.fastp.log"
        )

    # 3) minimap2 -> BAM (sorted+indexed)
    aln_bam = bc_out / f"{bc}.aln.bam"
    if not aln_bam.exists():
        cmd = (
            f"minimap2 -ax map-ont -t {threads} {shlex.quote(str(ref_fa))} {shlex.quote(str(filt_fq))} "
            f"| samtools view -bS -F 4 -@ {threads} - "
            f"| samtools sort -@ {threads} -o {shlex.quote(str(aln_bam))} -"
        )
        run(cmd, log_file=log_dir / f"{bc}.align.log")
        run(f"samtools index {shlex.quote(str(aln_bam))}")

    # 4) optional region filter
    region_bam = bc_out / f"{bc}.region.bam"
    if regions_file:
        if not region_bam.exists():
            run(
                f"samtools view -@ 8 -b -L {shlex.quote(str(regions_file))} "
                f"-o {shlex.quote(str(region_bam))} {shlex.quote(str(aln_bam))}"
            )
            run(f"samtools index {shlex.quote(str(region_bam))}")
    else:
        region_bam = aln_bam  # pass-through

    # 5) optional iVar trim
    if do_ivar:
        if not primer_bed:
            raise ValueError("--primer-bed required when --do-ivar-trim is set")
        ivar_prefix = bc_out / f"{bc}.primerclipped"
        ivar_bam = bc_out / f"{bc}.primerclipped.bam"
        if not ivar_bam.exists():
            run(
                " ".join([
                    "ivar trim",
                    f"-i {shlex.quote(str(region_bam))}",
                    f"-b {shlex.quote(str(primer_bed))}",
                    "-e -m 1 -s 4 -q 0",
                    f"-p {shlex.quote(str(ivar_prefix))}"
                ]),
                log_file=log_dir / f"{bc}.ivar_trim.log"
            )
        sort_input = ivar_bam
    else:
        sort_input = region_bam

    # 6) sort/index post-trim (or region)
    sorted_bam = bc_out / f"{bc}.sorted.bam"
    if not sorted_bam.exists():
        run(f"samtools sort -@ {threads} -o {shlex.quote(str(sorted_bam))} {shlex.quote(str(sort_input))}")
        run(f"samtools index {shlex.quote(str(sorted_bam))}")

    # 7) per-base cap (<= cap_depth) or pass
    calls_bam = bc_out / f"{bc}.calls.bam"
    if cap_depth and cap_depth > 0:
        if not calls_bam.exists():
            cap_bam_per_base(sorted_bam, calls_bam, cap_depth, regions_bed=regions_file, seed=cap_seed)
    else:
        # link/copy to keep naming consistent
        if not calls_bam.exists():
            run(f"ln -f {shlex.quote(str(sorted_bam))} {shlex.quote(str(calls_bam))}", check=False)
            run(f"ln -f {shlex.quote(str(sorted_bam))}.bai {shlex.quote(str(calls_bam))}.bai", check=False)

    # 8) coverage
    cov_txt = bc_out / f"{bc}.coverage.txt"
    if not cov_txt.exists():
        run(f"samtools coverage {shlex.quote(str(calls_bam))} > {shlex.quote(str(cov_txt))}")

    # 9) Clair3
    clair_out = bc_out / "clair3"
    clair_out.mkdir(exist_ok=True)
    merged_vcf = clair_out / "merge_output.vcf.gz"
    if not merged_vcf.exists():
        run(
            " ".join([
                "run_clair3.sh",
                f"-b {shlex.quote(str(calls_bam))}",
                f"-f {shlex.quote(str(ref_fa))}",
                f"-m {shlex.quote(str(clair3_model_dir))}",
                f"-p {shlex.quote(clair3_preset)}",
                f"-t {threads}",
                f"-o {shlex.quote(str(clair_out))}",
                "--include_all_ctgs",
            ]),
            log_file=log_dir / f"{bc}.clair3.log"
        )
        if not merged_vcf.exists():
            raise RuntimeError(f"Clair3 did not produce {merged_vcf}")

    # 10) snpEff
    ann_vcf = bc_out / f"{bc}.ann.vcf"
    if not ann_vcf.exists():
        run(f"snpEff -Xmx24G -v {shlex.quote(snpeff_db)} {shlex.quote(str(merged_vcf))} > {shlex.quote(str(ann_vcf))}")

    print(f"[OK] {bc} done.")

# ---------- CLI ----------

def parse_args():
    p = argparse.ArgumentParser(description="Barcode-aware ONT Pf panel pipeline (no Snakemake).")
    p.add_argument("--input-dir", required=True, type=Path, help="Folder containing barcode*/ subdirs")
    p.add_argument("--output-dir", required=True, type=Path, help="Where per-barcode outputs go")
    p.add_argument("--reference-fasta", required=True, type=Path, help="Reference FASTA (3D7)")
    p.add_argument("--regions-file", type=Path, default=None, help="BED/intervals for samtools -L (optional)")
    p.add_argument("--do-ivar-trim", action="store_true", help="Enable iVar primer trimming")
    p.add_argument("--primer-bed", type=Path, default=None, help="iVar BED (required if --do-ivar-trim)")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--fastp-min-qual", type=int, default=10)
    p.add_argument("--fastp-min-len", type=int, default=1000)
    p.add_argument("--fastp-max-len", type=int, default=2500)
    p.add_argument("--cap-depth", type=int, default=0, help="Per-base cap; 0 disables")
    p.add_argument("--cap-seed", type=int, default=42)
    p.add_argument("--clair3-model-dir", required=True, type=Path)
    p.add_argument("--clair3-preset", default="ont", choices=["ont", "hifi", "clr"])
    p.add_argument("--snpeff-db", required=True, help="e.g., pf3d7-61")
    p.add_argument("--barcode-prefix", default="barcode", help="Folder name prefix to detect barcodes")
    return p.parse_args()

def main():
    args = parse_args()
    # Preconditions
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ensure_fai(args.reference_fasta)

    barcodes = discover_barcodes(args.input_dir, args.barcode_prefix)
    if not barcodes:
        print(f"[ERROR] No '{args.barcode_prefix}*' dirs with FASTQs in {args.input_dir}", file=sys.stderr)
        sys.exit(2)
    print(f"[INFO] Found {len(barcodes)} barcodes")

    for bc_dir in barcodes:
        try:
            process_barcode(
                bc_dir=bc_dir,
                out_root=args.output_dir,
                threads=args.threads,
                ref_fa=args.reference_fasta,
                regions_file=args.regions_file,
                do_ivar=args.do_ivar_trim,
                primer_bed=args.primer_bed,
                fastp_q=args.fastp_min_qual,
                fastp_minlen=args.fastp_min_len,
                fastp_maxlen=args.fastp_max_len,
                cap_depth=args.cap_depth,
                cap_seed=args.cap_seed,
                clair3_model_dir=args.clair3_model_dir,
                clair3_preset=args.clair3_preset,
                snpeff_db=args.snpeff_db,
            )
        except Exception as e:
            print(f"[ERROR] {bc_dir.name}: {e}", file=sys.stderr)
            # continue to next barcode rather than aborting all
            continue

    print("[DONE] All barcodes processed.")

if __name__ == "__main__":
    main()

# EvidentialGenePipe

A config-driven pipeline that turns PacBio **FLNC IsoSeq reads** + a **reference genome**
+ a **reference annotation** into a **non-redundant coding transcriptome** with
per-transcript annotation of *reference* vs *novel* genes.

The final GFF3 is in **transcript coordinates** (seqid = transcript ID), not genome
coordinates. A genome-coordinate GFF is written alongside it purely as a cross-check.

Everything is CDS-focused: EvidentialGene classifies on coding sequence, and the
reference comparison is done on the transcript's ORF interval only.

---

## 1. Install

```bash
conda env create -f environment.yml
conda activate evigene_pipe
bash install_evigene.sh                 # downloads evigene23jul15.tar into ./software
```

The bioconda `isoseq` package (for `cluster2`) is linux-64 only, which is fine for RDAC
but will fail to solve on macOS - drop it from `environment.yml` and set
`tracks.cluster.enabled: false` if you ever need to run this on a Mac.

`install_evigene.sh` prints the `evigene_home` path and checks every external tool.
EvidentialGene is not on conda; `fastanrdb` (which tr2aacds needs) comes from the
bioconda `exonerate` package.

## 2. Configure

```bash
cp config/config.example.yaml config/Em2.yaml
$EDITOR config/Em2.yaml            # paths, threads, thresholds - nothing is hardcoded
python -m egpipe check --config config/Em2.yaml
```

`check` validates the config, confirms every executable is on PATH, prints the exact
`tr2aacds` command line that will run, and shows which stages are already done.

## 3. Run

Locally / on an interactive node:

```bash
python -m egpipe run --config config/Em2.yaml
```

Resume / rerun single stages:

```bash
python -m egpipe run --config config/Em2.yaml --from 04          # stage 04 onwards
python -m egpipe run --config config/Em2.yaml --only 05          # just reclassify
python -m egpipe run --config config/Em2.yaml --from 02 --force  # ignore checkpoints
```

Stages 03, 05 and 06 are cheap and always recomputed, so changing a `classify:` or
`transcriptome:` threshold and rerunning `--from 03` is fast. Stage 04 automatically
reruns if the transcript fasta is newer than the alignment.

---

## Stages

| # | What it does | Main outputs |
|---|---|---|
| 01 | `isoseq cluster2` on the FLNC BAM, IDs renamed, cluster size (`is:i`) kept as read support; reference CDS/protein extracted with `gffread` | `01_inputs/<SAMPLE>.fa`, `clustered.bam`, `flnc.idmap.tsv`, `ref_cds.raw.fa` |
| 02 | `tr2aacds4.pl` on the merged input | `02_evigene/okayset/`, `*.trclass`, `*.pubids` |
| 03 | Harvest okayset into a clean transcriptome; pull ORF offsets out of the evigene headers; store every transcript ORF-sense | `03_transcriptome/*.transcripts.fa`, `*.orf_table.tsv` |
| 04 | `minimap2 -x splice:hq -u f` transcripts -> genome | `04_mapping/*.genome.bam` |
| 05 | Classify reference vs novel; count FLNC support | `05_classify/*.annotation_table.tsv` |
| 06 | Write the transcript-coordinate GFF3 and the final deliverables | `06_final/` |

### Final outputs (`06_final/`)

* `<SAMPLE>.transcriptome.fa` - the unique transcript set
* `<SAMPLE>.cds.fa`, `<SAMPLE>.protein.fa`
* `<SAMPLE>.transcriptome.gff3` - **transcript coordinates**: one `gene`/`mRNA`/`exon`
  per transcript spanning `1..length`, the `CDS` at the evigene ORF offsets, and
  5'/3' UTR features. All annotation lives in the attributes (`class`, `ref_gene`,
  `ref_evidence`, `evigene_locus`, `aa_len`, `completeness`, `n_flnc_reads`,
  `genome_locus`, `genome_identity`).
* `<SAMPLE>.genome.gff3` - the same transcripts projected back onto the assembly
* `<SAMPLE>.annotation_table.tsv` - one row per transcript, all evidence columns
* `<SAMPLE>.summary.txt`, `<SAMPLE>.reference_genes_missing.txt`

### Classification

| class | meaning |
|---|---|
| `reference` | the transcript's CDS blocks overlap a reference coding gene's CDS exons on the same strand (`ref_evidence=cds_overlap`), or the transcript shares an evigene locus with a reference-CDS sequence (`evigene_locus` / `reference_cds_sequence`) |
| `novel_intergenic` | maps confidently, no reference coding gene at that locus |
| `novel_intronic` | maps inside an annotated coding gene's span but touches none of its CDS exons |
| `novel_antisense` | CDS blocks overlap reference CDS exons on the opposite strand (`antisense_to_gene` names the gene) |
| `novel_unmapped` | no alignment, or below `min_query_cov` / `min_identity` |

Reference genes that no transcript recovered are appended to the final set as
CDS-only entries with `ref_evidence=reference_only_no_flnc`, so the transcriptome
covers the whole reference gene set even where this sample's IsoSeq did not.

Direct genomic evidence wins: a transcript with antisense overlap is never promoted to
`reference` by locus sharing.

---

## Things worth knowing before you trust the numbers

* **`isoseq cluster2` is what keeps stage 02 tractable.** tr2aacds ends in an all-vs-all
  self-`blastn` of the surviving CDS, so the input size is the wall-clock and memory
  driver. cluster2 collapses reads into transcript consensuses first, and `--singletons`
  keeps the reads that formed no cluster, so nothing is thrown away. On the synthetic
  test set 35 reads become 15 transcripts; on a real sample expect a far larger factor.
  Set `tracks.cluster.enabled: false` to feed raw reads instead (identical reads are
  still collapsed, but that is a much weaker reduction). If stage 02 still exceeds the
  83 GB `defq` ceiling, re-render the sbatch scripts with `--partition bigmem`.
* **Chunking is the scaling path, and it is not free unless you carry drops forward.**
  A deep Kinnex run can leave ~1M clustered transcripts and the self-blastn grows as
  ~N^2. evigene documents running tr2aacds on subsets then re-running on the okay
  outputs, and says its classifier is deterministic with respect to how inputs are
  grouped. That claim holds for a fixed input set but *not* across chunking, because a
  chunk **drops** transcripts before any later stage sees them. Measured on the test set
  (32 raw reads, 4 chunks, single-pass = 30 transcripts):

  | mode | 1 level | 2 levels | 4 levels |
  |---|---|---|---|
  | `chunk_carry_forward: okay` | 29 (antisense 2/3) | 27 (antisense **0**/3) | 27 (antisense **0**/3) |
  | `chunk_carry_forward: nonredundant` | **30, identical** | **30, identical** | **30, identical** |

  The losses are not random: they concentrate in transcripts evigene drops as `parthi`
  or other high-identity alternate calls, which is exactly where antisense transcripts
  live. Carrying those forward fixes it. Only `perfectdup` and `perfectfrag` - sequences
  byte-identical to, or fully contained in, another sequence in the same chunk - are
  discarded inside a chunk, because nothing can be lost by doing so. Everything else is
  decided by the final combine, with the whole set in view.

  `evigene.chunks: auto` splits into `02_evigene/levelN/chunkNN/`, repeats until the
  surviving set fits `chunk_target`, then combines in `02_evigene/final/`. Sequences are
  renamed by the pipeline (`L1c1t0000001`, ...) with the mapping in
  `02_evigene/oid_chain.tsv`, because evigene truncates `-species` to 6 characters when
  building public IDs, so per-chunk species names collide and the merge aborts on
  duplicate IDs. `oids` in the final table still resolve to original cluster IDs.
* **Deciding what to do about singletons.** `--singletons` keeps every read that formed
  no cluster, which is where rare and novel transcripts live but also where most of the
  bulk and most of the error is. The cluster size is in `01_inputs/flnc.idmap.tsv`
  (`dupcount` column), so the singleton fraction is one `awk` away before you commit to
  a long run. Dropping them is a sensitivity decision, not a performance one - chunking
  solves the performance side without discarding evidence.
* **cluster2 hangs on a BAM without a PacBio read group.** A uBAM built by
  `samtools import` (or anything else that drops the `@RG` carrying `PL:PACBIO` and
  `DS:READTYPE=CCS`) makes `isoseq cluster2` sit forever with no output and no error.
  Stage 01 checks the header and fails with a message instead. FLNC BAMs as delivered
  by SMRT Link are fine.
* **Clustering costs you evigene's replication evidence.** `-aconsensus` treats the same
  CDS appearing in several input sequences as support. Once cluster2 has merged those
  reads into one consensus that signal is gone, which makes evigene more likely to drop
  ORFs under ~100 aa that also lack homology evidence (`evgclass=smallaa` in the
  dropset). Two knobs: `evigene.small_aa_max` (`-SMALLAAMAX`) lowers the threshold and
  `evigene.keep_small_orfs` (`-smallclass`) keeps small ORFs outright. Supplying
  homology via `-ablastab` in `evigene.extra_args` is the principled fix if small gene
  families matter to you.
* **Sequencing error still inflates the alternate count**, though far less after
  clustering. The locus count (`evigene_loci_total`) is the meaningful number, not the
  transcript count. Set `transcriptome.include_alternates: false` for one representative
  per locus, or raise `evigene.p_heterozygosity` (2-3) to collapse more aggressively.
* **Strand handling.** minimap2's own `-u` default is `n` (ignore GT-AG), which finds
  measurably worse intron boundaries - on the test set only 2 of 13 alignments
  reproduced the true exon structure under `-u n` versus 9 under `-u f` or `-u b`. The
  pipeline sets `-u f`, the documented Iso-Seq setting, which is only valid if the
  stored transcript is sense-oriented. Stage 03 guarantees that: evigene usually
  reverse-complements a minus-strand ORF itself (`revc=1`), but it can also report a
  reversed `offs=` while leaving the sequence alone, so stage 03 flips those and logs
  how many. Flipping cannot change a classification - reversing the sequence flips both
  the alignment strand and the ORF strand, and the genomic strand is their XOR.
  `mapping.strand_mode: "b"` remains available for mixed-orientation input. The cost of
  getting this wrong is measurable: feeding reverse-complemented transcripts to `-uf`
  drops exact exon-structure recovery from 9/13 to 5/13 on the test set, which is what
  the stage-03 reorientation guarantee prevents. PacBio's own ISOSEQ preset hardcodes
  `-uf`, which is independent corroboration of the default.
* **Why minimap2 directly and not pbmm2.** pbmm2 is PacBio's own minimap2 wrapper, so it
  is the obvious thing to reach for. Its `ISOSEQ` preset expands to
  `-k 15 -w 5 -u -o 6 -O 24 -e 1 -E 0 -A 1 -B 4 -z 200 -Z 100 -r 200000 -g 2000 -C 5
  -G 200000`, which is `minimap2 -x splice:hq` spelled out, and it also applies `-uf`.
  Verified rather than assumed: on the test transcripts pbmm2 ISOSEQ and
  `minimap2 -x splice:hq -uf` produce byte-identical alignments (13/13 records agreeing
  on contig, position, strand and every exon block), and the same holds on deliberately
  reverse-complemented input, which is the case where `-u` actually decides anything.
  So there is no accuracy argument either way. Two reasons to stay with minimap2 here:

  1. pbmm2's real value - carrying PacBio BAM tags, read groups and the `.pbi` through
     the alignment - applies to *reads*. Stage 04's input is evigene's derived transcript
     fasta, which has no PacBio metadata left; pbmm2 itself warns "Input is FASTA".
  2. pbmm2 filters by default (`-y 70` minimum gap-compressed identity, `-l 50` minimum
     mapped length). Stage 05 wants the alignment metrics for *every* transcript so it
     can apply `classify.min_identity` / `min_query_cov` itself and record the numbers in
     the annotation table. Filtering upstream would make a marginal alignment
     indistinguishable from no alignment at all, both landing as `novel_unmapped`.

  pbmm2 remains the right tool for aligning the FLNC reads themselves (as in a
  pbmm2 + `isoseq collapse` pipeline). Note also that `-u` means different things in the
  two programs: GT-AG strand search in minimap2, `--no-kmer-compression` in pbmm2.
* **`novel_unmapped` is not the same as novel.** Because the IsoSeq libraries come from a
  different subpopulation than the assembly, an unmapped transcript can be genuine
  novelty, population divergence, or an assembly gap. Treat that class as a follow-up
  list, not a result.
* **`n_flnc_reads` is a support count, not an expression estimate.** It comes from
  remapping the (deduplicated) FLNC reads to the final transcriptome with
  `minimap2 -x map-hifi --secondary=no`, counting forward-strand primary hits and
  weighting by how many identical reads each represents.
* **The reference CDS track is off by default, and that is deliberate.** Feeding the
  annotated CDS into tr2aacds alongside the reads looks appealing, but after cluster2
  the FLNC consensus CDS is frequently byte-identical to the annotated CDS. evigene's
  perfect-duplicate step then has to pick one, and it prefers the CDS-only sequence
  (100% coding beats ~70% coding) - so the full-length transcript is dropped as
  `perfectdup` and the reference gene loses its UTRs. On the test set that silently
  removed refgene1 and refgene4's real transcripts. Reference genes are instead added
  after classification via `transcriptome.add_unrecovered_reference`, which only fills
  the genuine gaps. `tracks.include_reference_cds: true` still works if you want the
  old behaviour.
* **Sample naming.** evigene derives its file names and public IDs from the input
  basename and breaks if it contains dots, so the pipeline sanitises `sample` to
  alphanumerics for the evigene run (`Em2` -> `Em2EVm000001t1`).
* tr2aacds can exit 0 after an internal `FATAL`; stage 02 scans the evigene logs and
  fails loudly if the okayset is unusable.

---

## Finding paralogs the reference does not resolve

`tools/find_collapsed_paralogs.py` looks for transcriptome loci that share one reference
gene or one genomic region - candidates for a gene family the annotation represents with
a single model, a tandem duplicate annotated as one gene, or two paralogs collapsed by
the assembly.

    python tools/find_collapsed_paralogs.py --config config/Em2.yaml
    python tools/find_collapsed_paralogs.py --config config/Em2.yaml --min-reads 5 --include-novel

It runs after stage 06 and writes `07_paralogs/<SAMPLE>.paralog_groups.tsv` (one row per
group, ranked) and `.paralog_pairs.tsv` (the pairwise evidence).

**It compares transcripts, not loci, and that is deliberate.** evigene builds loci from
high-identity local CDS alignments, and a paralog at ~93% CDS identity still shares
plenty of >=98% stretches. On the test set a deliberately planted paralog was filed
*inside* its parent's locus as an `alt` rather than given a locus of its own, so
"two loci on one reference gene" would have missed the very case being looked for. The
`same_evigene_locus` column tells you which situation you are in.

The discriminator is substitutions versus exon structure. Alternate isoforms differ by
which exons they carry, so where they align they are near-identical - identity stays
~100% while coverage drops. Paralogs differ by substitutions along their length -
coverage stays high while identity falls:

| verdict | meaning |
|---|---|
| `collapsed_paralog` | high coverage, identity below `--same-gene-identity` - two genes on one reference locus |
| `isoform_or_fragment` | identity at or above the cutoff - one gene split into evigene alternates |
| `tandem_or_split` | CDS footprints barely overlap - distinct copies under one model, or one model over two genes |
| `partial_overlap` | protein alignment covers too little of the shorter sequence to judge |

On the synthetic test (a planted paralog at 93% CDS identity plus nine error-driven
alternates) the separation was clean: the paralog pair at 91.6% protein identity, every
alternate pair at 98.1-99.2%. That margin is narrower than it looks, because the test
reads carry 0.2% error; real cluster2 consensus alternates sit closer to 100%. Either way the script prints the identity distribution of every pair it compared to
stdout on each run - 1% bins, bars log-scaled (the isoform mode is usually orders of
magnitude taller than any paralog mode), empty bins kept so the trough between modes is
visible, with the current cutoff marked and quantiles underneath:

    protein identity of compared pairs (n=10, --same-gene-identity=98.0, bars log-scaled)
       91.0-92.0  ############                             1
       92.0-93.0                                           0
       ...
       98.0-99.0  ######################################## 8  <-- current cutoff
       99.0-100.0 ############                             1
      min 91.56   p1 91.56   p5 91.56   p50 98.81   max 99.20

Set `--same-gene-identity` from that trough rather than trusting the default of 98. A
unimodal distribution pressed against 100 means there are no strong paralog candidates
to find.

`07_paralogs/` holds four files (blast intermediates go in `07_paralogs/blast/`):

| file | what it is |
|---|---|
| `<SAMPLE>.paralog_groups.tsv` | one row per group, ranked by verdict then read support |
| `<SAMPLE>.paralog_pairs.tsv` | the pairwise evidence - this is the raw data for re-binning |
| `<SAMPLE>.identity_histogram.txt` | the printed histogram, plus the parameters the run used |
| `<SAMPLE>.identity_histogram.tsv` | `bin_low`, `bin_high`, `count` - the bars as numbers |

Missing numeric values are written as `NA` rather than empty, so `read.delim()` types
those columns numeric instead of character. In R:

```r
pairs <- read.delim("07_paralogs/Em2.paralog_pairs.tsv")

# re-bin at whatever resolution you want
hist(pairs$protein_identity, breaks = seq(50, 100, by = 0.25),
     xlab = "protein identity (%)", main = "compared pairs")

# the discriminator is two-dimensional, so the scatter is more informative
# than either margin: isoforms sit top-right, paralogs drop down the y axis
plot(pairs$protein_aln_cov, pairs$protein_identity,
     col = factor(pairs$verdict), pch = 20,
     xlab = "protein alignment coverage", ylab = "protein identity (%)")

# within-locus pairs are the ones evigene already called alternates
table(pairs$same_evigene_locus, pairs$verdict)
```

Caveats worth keeping in mind: the verdict is a heuristic over two numbers, not a
phylogenetic test; `--min-reads` defaults to 2 because a pair of singletons differing by
5% is far more likely two error-laden reads than two genes; and each group is capped at
`--max-per-group` transcripts by read support, so a locus with hundreds of alternates
does not dominate the pairwise blastp.

## Troubleshooting

**cluster2 is OOM-killed (`oom_kill event`, `rc=-9`) while `sacct` reports a MaxRSS of
well under a gigabyte.** These look contradictory and are not. cluster2 spools a
temporary BAM of everything it writes and takes the directory from `TMPDIR`, falling
back to `/tmp`. On a cluster whose nodes have no local scratch (RDAC declares
`TmpDisk=0`), `/tmp` is RAM-backed - and tmpfs pages are charged to the job's memory
cgroup while never appearing in the process's RSS. So the job is killed for memory while
`MaxRSS` shows the process itself using almost none, and raising `--mem` only moves the
ceiling because the temp file keeps growing with the library. On RDAC this is literal -
`findmnt -T /tmp` on a defq node reports `rootfs 45.5G` against 92G of RAM, i.e. the
node is diskless and its entire root lives in memory. The same cause can also
surface as a write failure instead of a kill (see the next entry). Stage 01 now sets
`TMPDIR` and refuses to start if it points at a tmpfs:

    scratch: /home/.../Em5/tmp (ext4 on /home, 52000.0 GB free)      # good
    ERROR  resources.tmpdir (...) is on tmpfs, which is RAM.         # refused

Confirm what your nodes do with `srun -p defq --pty findmnt -no FSTYPE,SIZE /tmp`.

**`[pbbam] BAM writer ERROR: could not write record / file: /tmp/tmp.NNNN.bam.tmp /
reason: Resource temporarily unavailable`.** cluster2 spools a temporary BAM while it
writes, and it takes the directory from `TMPDIR`, falling back to `/tmp`. On a cluster
whose compute nodes have no local scratch (RDAC declares `TmpDisk=0`), that fallback
puts a temp file holding every record of a 20M-read library on the node's `/tmp`, and it
dies part-way through - note the RSS stays small, so this looks nothing like an OOM even
though it is often mistaken for one. The pipeline now sets `TMPDIR`/`TMP`/`TEMP` to
`resources.tmpdir` (default `<outdir>/tmp`) for cluster2, gffread and the stage-04
samtools sort, and logs the path it used:

    grep "cluster2 scratch" <outdir>/<SAMPLE>.egpipe.log

Point `resources.tmpdir` at a filesystem with room for a copy of the library.

**`isoseq cluster2` is killed for exceeding memory, and lowering `--cpus-per-task`
changes nothing.** Almost certainly because the lower CPU count never reached the tool.
`resources.threads` in the config drives `-j`, so an sbatch asking for 4 CPUs still ran
cluster2 with `-j 32` if that is what the config said, and `--sort-threads` defaults to
`-j` ("open file handles per BAM"). The pipeline now caps `resources.threads` and
`resources.max_mem_mb` at `SLURM_CPUS_PER_TASK` and `SLURM_MEM_PER_NODE`/`_PER_CPU`, and
logs a warning when it does, so the sbatch is the single source of truth. Verify from
the job log - it prints the exact command:

    grep "RUN  isoseq cluster2" <outdir>/<SAMPLE>.egpipe.log

cluster2 exposes no memory cap, so the levers are, in order: match `-j` to the
allocation (now automatic); set `tracks.cluster.sort_threads` to 2-4; then move to a
bigger partition. Two numbers worth having before choosing:

    seff <jobid>                                   # MaxRSS of the failed job
    tail -20 <outdir>/01_inputs/cluster2.log       # which phase it died in

The stage-01 sbatch profile previously hardcoded 8 CPU / 32G, which was wishful for a
20M-read library; it now follows `resources.threads` and `resources.max_mem_mb`.

**`Can't locate FindBin.pm in @INC` (or another core module) from an evigene helper.**
tr2aacds runs its helpers as executables (`$EVIGENES/cdna_bestorf.pl -cdna ...`), never
as `perl <script>`, so the `#!` line picks the interpreter - not PATH, and not the conda
environment the pipeline was launched from. 265 of the stock scripts hardcode
`#!/usr/bin/perl`, which lands on the system perl; on RHEL/Rocky that perl is modular
and core modules such as `FindBin` ship in separate RPMs that are often not installed.
`install_evigene.sh` now rewrites those shebangs to `#!/usr/bin/env perl`, and both
`egpipe check` and stage 02 execute a helper up front so this fails in seconds instead
of partway into a multi-hour job. To fix an existing install by hand:

    find <evigene_home>/scripts -name '*.pl' -exec \
      sed -i '1s|^#!.*[/ ]perl.*$|#!/usr/bin/env perl|' {} +

Five `.pl` files have no valid perl shebang at all and are left alone; none of them is
in the tr2aacds call path.

**`evigene reported problems ... FATAL ERR: failed bestorf_cds= ... nrec= 0` after a run
that clearly worked.** evigene opens its log with `>>`, so `<name>.tr2aacds.log`
accumulates every attempt in that directory. That FATAL comes from `failIfEmpty` right
after `get_bestorf`, which is step 0 - a run reporting it died immediately and produced
nothing, so it cannot be the run that just spent hours building a pubids table. It is
almost always left over from an earlier failed attempt (the FindBin problem above
produces exactly this, because the helper writes no CDS). Stage 02 now records each
log's size before launching and only scans the bytes that run appends, so stale errors
are no longer attributed to a successful run. To check a log by hand, compare line
numbers:

    grep -nE "EvidentialGene tr2aacds|FATAL" <name>.tr2aacds.log

If the last FATAL sits above the last `EvidentialGene tr2aacds` banner, it predates the
successful run.

## Tests

```bash
python tests/make_test_data.py --outdir /tmp/testdata        # synthetic genome/GFF/FLNC BAM
# point a config at /tmp/testdata, then:
python -m egpipe run --config /tmp/test_config.yaml
python tests/verify_outputs.py --outdir /tmp/testrun --sample TestSp \
    --expected /tmp/testdata/expected_classes.tsv
```

`make_test_data.py` builds a genome with five annotated genes (one on the minus strand,
one with no FLNC support), plus an intergenic novel gene, a gene inside an intron, an
antisense transcript, and a transcript absent from the genome. The FLNC BAM it writes
carries a proper PacBio `@RG` so `isoseq cluster2` will run on it.

`test_reorient.py` checks the ORF-sense guarantee that `-u f` depends on.

`verify_outputs.py` re-extracts every CDS from the transcriptome fasta using the GFF
coordinates and compares it byte-for-byte with the evigene CDS fasta, checks start/stop
codons and internal stops, checks that CDS+UTR tile each transcript exactly once, and
reports the class assigned to each planted origin.

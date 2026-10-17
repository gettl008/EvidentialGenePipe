#!/usr/bin/env python3
"""Check that stage 03 stores every transcript ORF-sense.

evigene normally reverse-complements a transcript whose best ORF is on the minus
strand (revc=1 in the header), but it can also report a reversed `offs=` while
leaving the sequence alone. Stage 04 uses `minimap2 -u f`, which assumes the
stored sequence is sense, so stage 03 has to fix those up.

  python tests/test_reorient.py
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from egpipe.config import load_config          # noqa: E402
from egpipe.stages import s03_harvest          # noqa: E402
from egpipe.utils import iter_fasta, revcomp, setup_logging  # noqa: E402

STOPS = ("TAA", "TAG", "TGA")


def main(tmp: Path = Path("/tmp/egpipe_reorient_test")) -> int:
    setup_logging(None, verbose=False)
    shutil.rmtree(tmp, ignore_errors=True)
    (tmp / "evg" / "okayset").mkdir(parents=True)

    orf = "ATG" + "GCT" * 20 + "TAA"
    lead, trail = "AAAA", "CCCCCC"
    seq = lead + revcomp(orf) + trail
    lo, hi = len(lead) + 1, len(lead) + len(orf)
    assert revcomp(seq[lo - 1:hi]) == orf

    fa = tmp / "evg" / "okayset" / "T.okay.mrna"
    fa.write_text(
        f">Tx1 type=mRNA; aalen={len(orf)//3-1},70%,complete; clen={len(seq)}; "
        f"strand=-; offs={hi}-{lo};\n{seq}\n"
    )
    (tmp / "c.yaml").write_text(
        f"sample: T\n"
        f"inputs: {{flnc_fasta: [{fa}], genome: {fa}, annotation: {fa}}}\n"
        f"outdir: {tmp}/out\n"
        f"evigene: {{home: {tmp}}}\n"
        f"tracks: {{cluster: {{enabled: false}}}}\n"
    )
    # config validation wants the evigene script to exist
    (tmp / "scripts" / "prot").mkdir(parents=True, exist_ok=True)
    (tmp / "scripts" / "prot" / "tr2aacds4.pl").write_text("#\n")

    cfg = load_config(tmp / "c.yaml")
    res = s03_harvest.run_stage(cfg, tmp / "evg")
    rec = res["records"]["Tx1"]
    out_seq = {n: s for n, _d, s in iter_fasta(res["transcripts_fa"])}["Tx1"]
    sub = out_seq[rec["cds_start"] - 1:rec["cds_end"]]

    checks = {
        "sequence reverse-complemented": out_seq == revcomp(seq),
        "cds_strand is +": rec["cds_strand"] == "+",
        "CDS slice equals the original ORF": sub == orf,
        "starts with ATG": sub.startswith("ATG"),
        "ends in a stop codon": sub[-3:] in STOPS,
        "length divisible by 3": len(sub) % 3 == 0,
        "no internal stop": all(sub[i:i + 3] not in STOPS
                                for i in range(0, len(sub) - 3, 3)),
    }
    for name, ok in checks.items():
        print(f"  {'ok  ' if ok else 'FAIL'} {name}")
    print(f"\n  CDS at {rec['cds_start']}-{rec['cds_end']} of {rec['length']} nt")
    failed = [k for k, v in checks.items() if not v]
    print("\nRESULT:", "PASS" if not failed else f"FAIL ({len(failed)})")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

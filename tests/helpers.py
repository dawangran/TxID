from __future__ import annotations

from pathlib import Path

from txid.workflow import initialize_registry

REFERENCE_GTF = """\
chr1\tref\texon\t100\t199\t.\t+\t.\tgene_id "g1"; transcript_id "tx_ref1"; gene_name "GENE1";
chr1\tref\texon\t300\t399\t.\t+\t.\tgene_id "g1"; transcript_id "tx_ref1"; gene_name "GENE1";
chr1\tref\texon\t500\t599\t.\t+\t.\tgene_id "g2"; transcript_id "tx_ref2";
chr1\tref\texon\t700\t799\t.\t+\t.\tgene_id "g2"; transcript_id "tx_ref2";
chr1\tref\texon\t1100\t1199\t.\t-\t.\tgene_id "g3"; transcript_id "tx_ref3";
chr1\tref\texon\t900\t999\t.\t-\t.\tgene_id "g3"; transcript_id "tx_ref3";
chr1\tref\texon\t1300\t1349\t.\t+\t.\tgene_id "gA"; transcript_id "tx_A";
chr1\tref\texon\t1450\t1499\t.\t+\t.\tgene_id "gA"; transcript_id "tx_A";
chr1\tref\texon\t1400\t1440\t.\t+\t.\tgene_id "gB"; transcript_id "tx_B";
chr1\tref\texon\t1550\t1599\t.\t+\t.\tgene_id "gB"; transcript_id "tx_B";
chr1\tref\texon\t1700\t1799\t.\t+\t.\tgene_id "gSE"; transcript_id "tx_SE";
"""

CALLER_GTF = """\
chr1\tcaller\texon\t300\t399\t42\t+\t.\tgene_id "up_g1"; transcript_id "alt_known"; note "known\\\"quote";
chr1\tcaller\texon\t100\t199\t41\t+\t.\tgene_id "up_g1"; transcript_id "alt_known"; note "known\\\"quote";
chr1\tcaller\texon\t300\t410\t.\t+\t.\tgene_id "up_g1"; transcript_id "end_variant";
chr1\tcaller\texon\t110\t199\t.\t+\t.\tgene_id "up_g1"; transcript_id "end_variant";
chr1\tcaller\texon\t300\t399\t.\t+\t.\tgene_id "up_g1"; transcript_id "one_bp_splice";
chr1\tcaller\texon\t100\t200\t.\t+\t.\tgene_id "up_g1"; transcript_id "one_bp_splice";
chr2\tcaller\texon\t100\t150\t.\t+\t.\tgene_id "new_upstream"; transcript_id "new_chr2";
chr2\tcaller\texon\t250\t300\t.\t+\t.\tgene_id "new_upstream"; transcript_id "new_chr2";
chr1\tcaller\texon\t100\t399\t.\t-\t.\tgene_id "antisense_upstream"; transcript_id "antisense_single";
chr1\tcaller\texon\t1320\t1360\t.\t+\t.\tgene_id "bridge"; transcript_id "ambiguous_AB";
chr1\tcaller\texon\t1500\t1570\t.\t+\t.\tgene_id "bridge"; transcript_id "ambiguous_AB";
chr1\tcaller\texon\t510\t790\t.\t+\t.\tgene_id "retained"; transcript_id "retained_intron";
chr1\tcaller\texon\t150\t180\t.\t+\t.\tgene_id "readthrough"; transcript_id "readthrough_1";
chr1\tcaller\texon\t720\t750\t.\t+\t.\tgene_id "readthrough"; transcript_id "readthrough_1";
chr1\tcaller\texon\t900\t999\t.\t-\t.\tgene_id "up_g3"; transcript_id "alt_minus";
chr1\tcaller\texon\t1100\t1199\t.\t-\t.\tgene_id "up_g3"; transcript_id "alt_minus";
"""


def write_inputs(directory: Path) -> tuple[Path, Path, Path]:
    fasta = directory / "reference.fa"
    reference = directory / "reference.gtf"
    caller = directory / "caller.gtf"
    fasta.write_text(
        ">chr1 synthetic\n" + ("ACGT" * 500) + "\n>chr2\n" + ("TGCA" * 250) + "\n",
        encoding="utf-8",
    )
    reference.write_text(REFERENCE_GTF, encoding="utf-8")
    caller.write_text(CALLER_GTF, encoding="utf-8")
    return fasta, reference, caller


def make_registry(directory: Path) -> tuple[Path, Path]:
    fasta, reference, caller = write_inputs(directory)
    database = directory / "txid.sqlite"
    initialize_registry(
        database,
        fasta=fasta,
        assembly_name="synthetic-v1",
        annotation=reference,
        annotation_name="ref-v1",
    )
    return database, caller


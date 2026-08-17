version 1.1

workflow TxIDBatch {
  input {
    Array[File] gtfs
    File reference_gtf
    File reference_fasta
    String assembly
    String annotation_name = "reference"
    File? contig_aliases
    Int? fuzzy_splice_tolerance_bp
    Int? fuzzy_end_tolerance_bp

    String docker_image = "dawang02/txid:0.1.3-jupyter"
    Int cpu = 1
    Int memory_gb = 4
    Int disk_gb = 20
  }

  call RunTxIDBatch {
    input:
      gtfs = gtfs,
      reference_gtf = reference_gtf,
      reference_fasta = reference_fasta,
      assembly = assembly,
      annotation_name = annotation_name,
      contig_aliases = contig_aliases,
      fuzzy_splice_tolerance_bp = fuzzy_splice_tolerance_bp,
      fuzzy_end_tolerance_bp = fuzzy_end_tolerance_bp,
      docker_image = docker_image,
      cpu = cpu,
      memory_gb = memory_gb,
      disk_gb = disk_gb
  }

  output {
    File registry = RunTxIDBatch.registry
    File catalog = RunTxIDBatch.catalog
    File manifest = RunTxIDBatch.manifest
    File initialization_summary = RunTxIDBatch.initialization_summary
    File batch_summary = RunTxIDBatch.batch_summary
    File validation_summary = RunTxIDBatch.validation_summary
    Array[File] rewritten_gtfs = RunTxIDBatch.rewritten_gtfs
    Array[File] mappings = RunTxIDBatch.mappings
  }
}

task RunTxIDBatch {
  input {
    Array[File] gtfs
    File reference_gtf
    File reference_fasta
    String assembly
    String annotation_name
    File? contig_aliases
    Int? fuzzy_splice_tolerance_bp
    Int? fuzzy_end_tolerance_bp
    String docker_image
    Int cpu
    Int memory_gb
    Int disk_gb
  }

  File gtf_list = write_lines(gtfs)
  File annotation_name_file = write_lines([annotation_name])
  File assembly_file = write_lines([assembly])

  command <<<
    if [ "~{defined(fuzzy_splice_tolerance_bp)}" != "~{defined(fuzzy_end_tolerance_bp)}" ]; then
      echo "fuzzy mode requires both fuzzy tolerances" >&2
      exit 2
    fi

    python - \
      "~{gtf_list}" \
      "~{annotation_name_file}" \
      "~{assembly_file}" \
      import-manifest.tsv <<'PY'
    import csv
    import sys
    from pathlib import Path

    gtf_list, annotation_file, assembly_file, output = sys.argv[1:]

    def read_lines(path):
        return Path(path).read_text(encoding="utf-8").splitlines()

    gtfs = read_lines(gtf_list)
    annotation_names = read_lines(annotation_file)
    assemblies = read_lines(assembly_file)

    if not gtfs:
        raise SystemExit("gtfs must contain at least one input GTF")
    if len(annotation_names) != 1 or len(assemblies) != 1:
        raise SystemExit("annotation_name and assembly must each contain exactly one line")

    annotation_name = annotation_names[0]
    assembly = assemblies[0]
    for label, values in (
        ("annotation_name", [annotation_name]),
        ("assembly", [assembly]),
    ):
        for value in values:
            if not value or "\t" in value or "\n" in value or "\r" in value:
                raise SystemExit(
                    f"{label} values must be non-empty and cannot contain tabs or newlines"
                )

    with Path(output).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("input", "sample", "tool", "annotation_name", "format"))
        for index, gtf in enumerate(gtfs, start=1):
            filename = Path(gtf).name
            lowered = filename.lower()
            if lowered.endswith(".gtf.gz"):
                sample = filename[:-7]
            elif lowered.endswith(".gtf"):
                sample = filename[:-4]
            else:
                sample = Path(filename).stem
            if not sample:
                sample = f"input_{index:03d}"
            writer.writerow((gtf, sample, "unspecified", annotation_name, "gtf"))
    PY
    manifest_status=$?
    if [ "$manifest_status" -ne 0 ]; then
      exit "$manifest_status"
    fi

    annotation_name_value=$(cat "~{annotation_name_file}") || exit $?
    assembly_value=$(cat "~{assembly_file}") || exit $?

    if [ "~{defined(contig_aliases)}" = "true" ]; then
      txid init \
        --db txid.sqlite \
        --fasta "~{reference_fasta}" \
        --assembly "$assembly_value" \
        --annotation "~{reference_gtf}" \
        --annotation-name "$annotation_name_value" \
        --contig-aliases "~{contig_aliases}" \
        > txid.init.json || exit $?
    else
      txid init \
        --db txid.sqlite \
        --fasta "~{reference_fasta}" \
        --assembly "$assembly_value" \
        --annotation "~{reference_gtf}" \
        --annotation-name "$annotation_name_value" \
        > txid.init.json || exit $?
    fi

    mkdir outputs || exit $?
    if [ "~{defined(fuzzy_splice_tolerance_bp)}" = "true" ]; then
      txid batch \
        --db txid.sqlite \
        --manifest import-manifest.tsv \
        --output-dir outputs \
        --fuzzy-splice-tolerance "~{fuzzy_splice_tolerance_bp}" \
        --fuzzy-end-tolerance "~{fuzzy_end_tolerance_bp}" \
        > txid.batch.json || exit $?
    else
      txid batch \
        --db txid.sqlite \
        --manifest import-manifest.tsv \
        --output-dir outputs \
        > txid.batch.json || exit $?
    fi

    txid export --db txid.sqlite --catalog txid.catalog.tsv > txid.export.json || exit $?
    txid validate --db txid.sqlite > txid.validation.json || exit $?
  >>>

  output {
    File registry = "txid.sqlite"
    File catalog = "txid.catalog.tsv"
    File manifest = "import-manifest.tsv"
    File initialization_summary = "txid.init.json"
    File batch_summary = "txid.batch.json"
    File validation_summary = "txid.validation.json"
    Array[File] rewritten_gtfs = glob("outputs/*.txid.gtf")
    Array[File] mappings = glob("outputs/*.mapping.tsv")
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} SSD"
  }
}

version 1.1

workflow TxIDMultiAdd {
  input {
    File reference_fasta
    File reference_gtf
    Array[File] gtfs
    String assembly
    String annotation_name = "reference"
    String tool = "unspecified"
    File? contig_aliases
    Int? fuzzy_splice_tolerance_bp
    Int? fuzzy_end_tolerance_bp

    String docker_image = "dawang02/txid:0.1.3-jupyter"
    Int cpu = 1
    Int memory_gb = 4
    Int disk_gb = 20
  }

  call RunTxIDMultiAdd {
    input:
      reference_fasta = reference_fasta,
      reference_gtf = reference_gtf,
      gtfs = gtfs,
      assembly = assembly,
      annotation_name = annotation_name,
      tool = tool,
      contig_aliases = contig_aliases,
      fuzzy_splice_tolerance_bp = fuzzy_splice_tolerance_bp,
      fuzzy_end_tolerance_bp = fuzzy_end_tolerance_bp,
      docker_image = docker_image,
      cpu = cpu,
      memory_gb = memory_gb,
      disk_gb = disk_gb
  }

  output {
    File registry = RunTxIDMultiAdd.registry
    File catalog = RunTxIDMultiAdd.catalog
    File initialization_summary = RunTxIDMultiAdd.initialization_summary
    File multi_add_summary = RunTxIDMultiAdd.multi_add_summary
    File export_summary = RunTxIDMultiAdd.export_summary
    File validation_summary = RunTxIDMultiAdd.validation_summary
    Array[File] rewritten_gtfs = RunTxIDMultiAdd.rewritten_gtfs
    Array[File] mappings = RunTxIDMultiAdd.mappings
  }
}

task RunTxIDMultiAdd {
  input {
    File reference_fasta
    File reference_gtf
    Array[File] gtfs
    String assembly
    String annotation_name
    String tool
    File? contig_aliases
    Int? fuzzy_splice_tolerance_bp
    Int? fuzzy_end_tolerance_bp
    String docker_image
    Int cpu
    Int memory_gb
    Int disk_gb
  }

  File gtf_list = write_lines(gtfs)
  File assembly_file = write_lines([assembly])
  File annotation_name_file = write_lines([annotation_name])
  File tool_file = write_lines([tool])

  command <<<
    python - \
      "~{gtf_list}" \
      "~{reference_fasta}" \
      "~{reference_gtf}" \
      "~{assembly_file}" \
      "~{annotation_name_file}" \
      "~{tool_file}" \
      "~{defined(contig_aliases)}" \
      "~{contig_aliases}" \
      "~{defined(fuzzy_splice_tolerance_bp)}" \
      "~{fuzzy_splice_tolerance_bp}" \
      "~{defined(fuzzy_end_tolerance_bp)}" \
      "~{fuzzy_end_tolerance_bp}" <<'PY'
    import subprocess
    import sys
    from pathlib import Path

    (
        gtf_list,
        reference_fasta,
        reference_gtf,
        assembly_file,
        annotation_name_file,
        tool_file,
        aliases_defined,
        contig_aliases,
        fuzzy_splice_defined,
        fuzzy_splice_value,
        fuzzy_end_defined,
        fuzzy_end_value,
    ) = sys.argv[1:]

    def read_scalar(path, label):
        values = Path(path).read_text(encoding="utf-8").splitlines()
        if len(values) != 1 or not values[0] or "\t" in values[0]:
            raise SystemExit(
                f"{label} must be one non-empty line without tab characters"
            )
        return values[0]

    def optional_nonnegative(defined, value, label):
        if defined != "true":
            return None
        try:
            parsed = int(value)
        except ValueError as error:
            raise SystemExit(f"{label} must be an integer") from error
        if parsed < 0:
            raise SystemExit(f"{label} must be non-negative")
        return parsed

    def run_json(output, arguments):
        with Path(output).open("w", encoding="utf-8") as handle:
            subprocess.run(arguments, stdout=handle, check=True)

    gtfs = Path(gtf_list).read_text(encoding="utf-8").splitlines()
    if not gtfs:
        raise SystemExit("gtfs must contain at least one input GTF")
    if any(not Path(gtf).is_file() for gtf in gtfs):
        raise SystemExit("every localized gtfs entry must be a file")

    assembly = read_scalar(assembly_file, "assembly")
    annotation_name = read_scalar(annotation_name_file, "annotation_name")
    tool = read_scalar(tool_file, "tool")

    if fuzzy_splice_defined != fuzzy_end_defined:
        raise SystemExit("fuzzy mode requires both fuzzy tolerances")
    fuzzy_splice = optional_nonnegative(
        fuzzy_splice_defined,
        fuzzy_splice_value,
        "fuzzy_splice_tolerance_bp",
    )
    fuzzy_end = optional_nonnegative(
        fuzzy_end_defined,
        fuzzy_end_value,
        "fuzzy_end_tolerance_bp",
    )

    init_command = [
        "txid",
        "init",
        "--db",
        "txid.sqlite",
        "--fasta",
        reference_fasta,
        "--assembly",
        assembly,
        "--annotation",
        reference_gtf,
        "--annotation-name",
        annotation_name,
        "--format",
        "gtf",
    ]
    if aliases_defined == "true":
        init_command.extend(["--contig-aliases", contig_aliases])
    run_json("txid.init.json", init_command)

    Path("outputs").mkdir()
    multi_add_command = [
        "txid",
        "multi-add",
        "--db",
        "txid.sqlite",
        "--input",
        *gtfs,
        "--tool",
        tool,
        "--annotation-name",
        annotation_name,
        "--output-dir",
        "outputs",
        "--format",
        "gtf",
    ]
    if fuzzy_splice is not None and fuzzy_end is not None:
        multi_add_command.extend(
            [
                "--fuzzy-splice-tolerance",
                str(fuzzy_splice),
                "--fuzzy-end-tolerance",
                str(fuzzy_end),
            ]
        )
    run_json("txid.multi-add.json", multi_add_command)

    run_json(
        "txid.export.json",
        [
            "txid",
            "export",
            "--db",
            "txid.sqlite",
            "--catalog",
            "txid.catalog.tsv",
        ],
    )
    run_json(
        "txid.validation.json",
        ["txid", "validate", "--db", "txid.sqlite"],
    )
    PY
  >>>

  output {
    File registry = "txid.sqlite"
    File catalog = "txid.catalog.tsv"
    File initialization_summary = "txid.init.json"
    File multi_add_summary = "txid.multi-add.json"
    File export_summary = "txid.export.json"
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

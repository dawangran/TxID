PRAGMA foreign_keys = ON;

CREATE TABLE schema_meta (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE reference_context (
    id INTEGER PRIMARY KEY CHECK (id = 1),
    assembly_name TEXT NOT NULL,
    fingerprint TEXT NOT NULL UNIQUE,
    fasta_checksum TEXT NOT NULL,
    alias_checksum TEXT,
    created_utc TEXT NOT NULL
);

CREATE TABLE contig (
    name TEXT PRIMARY KEY,
    length INTEGER NOT NULL CHECK (length > 0),
    sequence_digest TEXT NOT NULL
);

CREATE TABLE contig_alias (
    alias TEXT PRIMARY KEY,
    primary_name TEXT NOT NULL REFERENCES contig(name)
);

CREATE TABLE structural_object (
    public_id TEXT PRIMARY KEY,
    family TEXT NOT NULL CHECK (family IN ('SC1', 'TF1', 'SE1')),
    public_digest TEXT NOT NULL,
    full_digest TEXT NOT NULL,
    canonical_json TEXT NOT NULL,
    UNIQUE (family, full_digest)
);

CREATE TABLE annotation_context (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    fingerprint TEXT NOT NULL,
    input_checksum TEXT NOT NULL,
    created_utc TEXT NOT NULL
);

CREATE TABLE reference_transcript (
    id INTEGER PRIMARY KEY,
    annotation_id INTEGER NOT NULL REFERENCES annotation_context(id),
    reference_gene_id TEXT,
    reference_transcript_id TEXT NOT NULL,
    contig TEXT NOT NULL REFERENCES contig(name),
    strand TEXT NOT NULL CHECK (strand IN ('+', '-')),
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    exons_json TEXT NOT NULL,
    attributes_json TEXT NOT NULL,
    splice_chain_id TEXT REFERENCES structural_object(public_id),
    form_id TEXT NOT NULL REFERENCES structural_object(public_id),
    UNIQUE (annotation_id, reference_transcript_id)
);

CREATE INDEX reference_transcript_structure
ON reference_transcript(annotation_id, form_id);

CREATE INDEX reference_transcript_location
ON reference_transcript(annotation_id, contig, strand, start, end);

CREATE TABLE import_manifest (
    id INTEGER PRIMARY KEY,
    annotation_id INTEGER NOT NULL REFERENCES annotation_context(id),
    sample TEXT NOT NULL,
    upstream_tool TEXT NOT NULL,
    source_label TEXT NOT NULL,
    input_checksum TEXT NOT NULL,
    options_json TEXT NOT NULL,
    software_version TEXT NOT NULL,
    created_utc TEXT NOT NULL,
    UNIQUE (annotation_id, sample, upstream_tool, input_checksum)
);

CREATE TABLE locus (
    accession INTEGER PRIMARY KEY AUTOINCREMENT,
    public_id TEXT NOT NULL UNIQUE,
    contig TEXT NOT NULL REFERENCES contig(name),
    strand TEXT NOT NULL CHECK (strand IN ('+', '-')),
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    status TEXT NOT NULL CHECK (status IN ('new_locus', 'ambiguous_gene')),
    gene_candidates_json TEXT NOT NULL,
    created_utc TEXT NOT NULL
);

CREATE INDEX locus_location ON locus(contig, strand, start, end);

CREATE TABLE observation (
    id INTEGER PRIMARY KEY,
    import_id INTEGER NOT NULL REFERENCES import_manifest(id),
    ordinal INTEGER NOT NULL,
    original_gene_id TEXT,
    original_transcript_id TEXT NOT NULL,
    output_gene_id TEXT NOT NULL,
    output_transcript_id TEXT NOT NULL,
    classification TEXT NOT NULL CHECK (
        classification IN ('known', 'novel_in_known_gene', 'ambiguous_gene', 'new_locus')
    ),
    gene_candidates_json TEXT NOT NULL,
    contig TEXT NOT NULL REFERENCES contig(name),
    strand TEXT NOT NULL CHECK (strand IN ('+', '-')),
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    exons_json TEXT NOT NULL,
    attributes_json TEXT NOT NULL,
    exon_attributes_json TEXT NOT NULL,
    source TEXT NOT NULL,
    source_format TEXT NOT NULL,
    exon_metadata_json TEXT NOT NULL,
    splice_chain_id TEXT REFERENCES structural_object(public_id),
    form_id TEXT NOT NULL REFERENCES structural_object(public_id),
    locus_accession INTEGER REFERENCES locus(accession),
    fuzzy_cluster_id TEXT,
    UNIQUE (import_id, original_transcript_id)
);

CREATE INDEX observation_form ON observation(form_id);
CREATE INDEX observation_output ON observation(output_transcript_id);

CREATE TABLE fuzzy_run (
    id INTEGER PRIMARY KEY,
    algorithm TEXT NOT NULL,
    splice_tolerance INTEGER NOT NULL CHECK (splice_tolerance >= 0),
    end_tolerance INTEGER NOT NULL CHECK (end_tolerance >= 0),
    input_fingerprint TEXT NOT NULL,
    software_version TEXT NOT NULL,
    created_utc TEXT NOT NULL,
    UNIQUE (algorithm, splice_tolerance, end_tolerance, input_fingerprint)
);

CREATE TABLE fuzzy_cluster (
    accession INTEGER PRIMARY KEY AUTOINCREMENT,
    public_id TEXT NOT NULL UNIQUE,
    fuzzy_run_id INTEGER NOT NULL REFERENCES fuzzy_run(id),
    signature TEXT NOT NULL,
    bridge_status TEXT NOT NULL CHECK (bridge_status IN ('unambiguous', 'ambiguous_bridge')),
    UNIQUE (fuzzy_run_id, signature)
);

CREATE TABLE fuzzy_member (
    fuzzy_cluster_accession INTEGER NOT NULL REFERENCES fuzzy_cluster(accession),
    form_id TEXT NOT NULL REFERENCES structural_object(public_id),
    PRIMARY KEY (fuzzy_cluster_accession, form_id)
);

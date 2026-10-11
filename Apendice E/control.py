"""
control.py

Integrated synonymous-recoding null test, matched-control selection,
control null test, and individual focus-control comparison for
Virus-Host codon-FCGR analyses.


WORKFLOW
--------

1. Load Virus and focus Host CDS sequences.

2. Construct a host_global synonymous-codon probability model
   from the supplied focus Host CDSs.

3. Run the synonymous-recoding null test for the focus genes.

4. Select matched background CDS controls using:
       - CDS length
       - GC content
       - GC3 content
       - amino-acid Jensen-Shannon distance

5. Construct a separate host_global synonymous-codon probability
   model from all selected control CDSs.

6. Run the synonymous-recoding null test for the selected controls.

7. Compare every selected control individually with its corresponding
   focus gene.

No averaging or median aggregation across matched controls is performed.


REQUIRED INPUTS
---------------

-v / --virus
    Virus FASTA file or directory containing Virus FASTA files.

-h / --host
    Focus Host FASTA file or directory containing focus Host FASTA files.

-t / --trna
    Plain-text tRNA copy-number table.

-j / --sij
    Plain-text file containing:
        s_GU
        s_UG
        s_IC
        s_IA

-i / --input
    Multi-FASTA background CDS file used for control selection.

-s / --seed
    Random seed.

-p / --perm
    Number of synonymous-recoding permutations.

-c / --controls
    Number of matched controls selected per focus gene.


NULL TEST
---------

Modes:
    raw
    penalized

Null model:
    host_global

Host-global synonymous pseudocount:
    1.0

Metrics:
    spearman_r
    pearson_r
    ssim
    competition_normalized_score

Standardized statistical implementations:
    Spearman:
        scipy.stats.spearmanr

    Pearson:
        scipy.stats.pearsonr

    SSIM:
        skimage.metrics.structural_similarity

    Benjamini-Hochberg:
        statsmodels.stats.multitest.multipletests(method="fdr_bh")


For amino acid a and synonymous codon c:

                         N(a,c) + 1
    P(c | a) = --------------------------------
               sum_j [N(a,c_j) + 1]


Focus and control synonymous probability models are estimated
independently.

Focus null model:
    Estimated from all supplied focus Host CDSs.

Control null model:
    Estimated from all selected control CDSs.


Within each null-test execution, one RNG is initialized before
the mode loop.

RAW mode consumes the first random values.

PENALIZED mode continues from the RNG state left after RAW.

The focus null test and control null test each initialize their
own RNG from the same supplied seed.


CONTROL SELECTION
-----------------

Minimum background CDS length:
    30 sense codons

Matching variables:
    1. CDS length
    2. GC content
    3. GC3 content
    4. amino-acid composition

Length distance:

    abs(log(control_codons / focus_codons))

GC distance:

    abs(control_GC - focus_GC)

GC3 distance:

    abs(control_GC3 - focus_GC3)

Amino-acid distance:

    Jensen-Shannon distance between amino-acid
    frequency distributions.

Jensen-Shannon zero handling:

    Exact zero-frequency categories contribute zero to their
    corresponding KL term.

    Boolean masks are used.

    No epsilon or pseudocount is added.

Each matching component is standardized independently using
population z-scores within the currently available candidate
pool for each focus gene.

Final matching score:

    z_length
    + z_GC
    + z_GC3
    + z_amino_acid
    + tiny seeded jitter

Lower scores indicate better matches.

Control selection uses its own RNG initialized from the supplied
seed and does not consume null-test random numbers.


FOCUS-CONTROL DIFFERENCES
-------------------------

All direct focus-control differences are absolute.

Observed:

    abs_observed_difference
        = abs(focus_observed - control_observed)

Null mean:

    abs_null_mean_difference
        = abs(focus_null_mean - control_null_mean)

Null median:

    abs_null_median_difference
        = abs(focus_null_median - control_null_median)


OUTPUTS
-------

The files are written directly beside this script.

null_test.xlsx

    Focus
        Null-test results for focus genes.

    Control
        Individual focus-control comparison results.

controls.xlsx

    Genes_control
        Matched-control selection report.
"""


# ============================================================
# IMPORTS
# ============================================================

import argparse
import math
import re
import time
import warnings

from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from scipy import stats
from statsmodels.stats.multitest import multipletests
from skimage.metrics import structural_similarity as ssim


# ============================================================
# GENERAL CONSTANTS
# ============================================================

DNA_BASES = "ACGT"

CODON_LENGTH = 3

MATRIX_SIZE = 8

NULL_WORKBOOK_NAME = "null_test.xlsx"

CONTROLS_WORKBOOK_NAME = "controls.xlsx"

CONTROL_SELECTION_SHEET = "Genes_control"

MIN_CONTROL_CODONS = 30

HOST_GLOBAL_PSEUDOCOUNT = 1.0

MODES = (
    "raw",
    "penalized",
)

METRICS = (
    "spearman_r",
    "pearson_r",
    "ssim",
    "competition_normalized_score",
)

FASTA_EXTENSIONS = {
    ".fa",
    ".fas",
    ".fasta",
    ".fna",
}


# ============================================================
# GENETIC CODE
# ============================================================

GENETIC_CODE = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",

    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",

    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",

    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",

    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",

    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",

    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",

    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",

    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",

    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",

    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",

    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",

    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",

    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",

    "TAC": "Y",
    "TAT": "Y",
    "TAA": "*",
    "TAG": "*",

    "TGC": "C",
    "TGT": "C",
    "TGA": "*",
    "TGG": "W",
}


ALL_CODONS = sorted(
    GENETIC_CODE
)


STOP_CODONS = {
    codon

    for codon, amino_acid
    in GENETIC_CODE.items()

    if amino_acid == "*"
}


SENSE_CODONS = sorted(
    codon

    for codon, amino_acid
    in GENETIC_CODE.items()

    if amino_acid != "*"
)


AA_ORDER = sorted(
    {
        amino_acid

        for amino_acid
        in GENETIC_CODE.values()

        if amino_acid != "*"
    }
)


SIJ_PARAMETERS = (
    "s_GU",
    "s_UG",
    "s_IC",
    "s_IA",
)


# ============================================================
# FCGR CORNERS
# ============================================================

CORNER_MAPPING = {
    "A": (0, 0),
    "C": (1, 0),
    "G": (0, 1),
    "T": (1, 1),
}


# ============================================================
# WOBBLE MODEL
# ============================================================

I34_CANDIDATE_ANTICODONS_RNA = {
    "AGC",
    "AAU",
    "AAG",
    "AGG",
    "ACG",
    "AGA",
    "AGU",
    "AAC",
}


WOBBLE_THIRD_BASES = {
    "C": ("G",),
    "A": ("U",),
    "G": ("C", "U"),
    "U": ("A", "G"),
    "I": ("U", "C", "A"),
}


PAIR_TO_SIJ_PARAMETER = {
    "G:U": "s_GU",
    "U:G": "s_UG",
    "I:C": "s_IC",
    "I:A": "s_IA",
}


# ============================================================
# FCGR COORDINATES
# ============================================================

def codon_coordinate(codon):

    x = 0
    y = 0

    for base in codon:

        dx, dy = (
            CORNER_MAPPING[
                base
            ]
        )

        x = (
            (x << 1)
            | dx
        )

        y = (
            (y << 1)
            | dy
        )

    return x, y


CODON_TO_COORDINATE = {
    codon:
        codon_coordinate(
            codon
        )

    for codon
    in ALL_CODONS
}


# ============================================================
# PATH HANDLING
# ============================================================

def script_directory():

    return (
        Path(__file__)
        .resolve()
        .parent
    )


def resolve_path(path_argument):

    path = Path(
        path_argument
    )

    if not path.is_absolute():

        path = (
            script_directory()
            / path
        )

    return (
        path.resolve()
    )


def validate_sequence_source(
    path,
    label
):

    if not path.exists():

        raise FileNotFoundError(
            f"{label} input not found:\n"
            f"{path}"
        )

    if not (
        path.is_file()
        or path.is_dir()
    ):

        raise ValueError(
            f"{label} input is neither "
            f"a file nor directory:\n"
            f"{path}"
        )


def validate_file(
    path,
    label
):

    if not path.exists():

        raise FileNotFoundError(
            f"{label} file not found:\n"
            f"{path}"
        )

    if not path.is_file():

        raise ValueError(
            f"{label} input must be a file:\n"
            f"{path}"
        )


# ============================================================
# FASTA DISCOVERY
# ============================================================

def find_fasta_files(folder):

    fasta_files = sorted(
        (
            path

            for path
            in folder.iterdir()

            if (
                path.is_file()
                and path.suffix.lower()
                in FASTA_EXTENSIONS
            )
        ),

        key=lambda path:
            path.name.lower()
    )

    if not fasta_files:

        raise ValueError(
            f"No FASTA files found in:\n"
            f"{folder}"
        )

    return fasta_files


# ============================================================
# FASTA READER
# ============================================================

def read_fasta_records(
    file_path
):

    records = []

    current_header = None

    sequence_parts = []

    with open(
        file_path,
        "r",
        encoding="utf-8-sig"
    ) as handle:

        for (
            line_number,
            line

        ) in enumerate(
            handle,
            start=1
        ):

            line = (
                line.strip()
            )

            if not line:
                continue

            if line.startswith(">"):

                if current_header is not None:

                    sequence = (
                        "".join(
                            sequence_parts
                        )
                        .upper()
                    )

                    if not sequence:

                        raise ValueError(
                            f"Empty FASTA sequence for "
                            f"'{current_header}' in:\n"
                            f"{file_path}"
                        )

                    records.append(
                        (
                            current_header,
                            sequence
                        )
                    )

                current_header = (
                    line[1:]
                    .strip()
                )

                if not current_header:

                    raise ValueError(
                        f"Empty FASTA header on line "
                        f"{line_number} in:\n"
                        f"{file_path}"
                    )

                sequence_parts = []

            else:

                if current_header is None:

                    raise ValueError(
                        f"Sequence data before first FASTA "
                        f"header on line {line_number} in:\n"
                        f"{file_path}"
                    )

                sequence_parts.append(
                    re.sub(
                        r"\s+",
                        "",
                        line
                    )
                )

    if current_header is not None:

        sequence = (
            "".join(
                sequence_parts
            )
            .upper()
        )

        if not sequence:

            raise ValueError(
                f"Empty FASTA sequence for "
                f"'{current_header}' in:\n"
                f"{file_path}"
            )

        records.append(
            (
                current_header,
                sequence
            )
        )

    if not records:

        raise ValueError(
            f"No FASTA records found in:\n"
            f"{file_path}"
        )

    return records


# ============================================================
# SEQUENCE SOURCE LOADER
# ============================================================

def load_sequence_source(
    source_path,
    role
):

    if source_path.is_file():

        fasta_files = [
            source_path
        ]

    else:

        fasta_files = (
            find_fasta_files(
                source_path
            )
        )

    sequences = []

    for file_path in fasta_files:

        records = (
            read_fasta_records(
                file_path
            )
        )

        multiple_records = (
            len(records) > 1
        )

        base_name = (
            file_path.stem
        )

        for (
            record_number,
            (
                header,
                sequence
            )

        ) in enumerate(
            records,
            start=1
        ):

            sequence_id = (
                header
                .split()[0]
            )

            if multiple_records:

                sequence_name = (
                    f"{base_name}_"
                    f"{record_number}"
                )

            else:

                sequence_name = (
                    base_name
                )

            sequences.append(
                {
                    "id":
                        sequence_id,

                    "name":
                        sequence_name,

                    "header":
                        header,

                    "sequence":
                        sequence,

                    "source_file":
                        file_path,

                    "role":
                        role,
                }
            )

    if not sequences:

        raise ValueError(
            f"No usable {role} FASTA records detected."
        )

    return sequences


# ============================================================
# CODON EXTRACTION FOR NULL TEST
# ============================================================

def extract_codons(sequence):

    sequence = (
        sequence.upper()
    )

    codons = []

    for position in range(
        0,
        len(sequence) - CODON_LENGTH + 1,
        CODON_LENGTH
    ):

        codon = sequence[
            position:
            position + CODON_LENGTH
        ]

        if (
            len(codon)
            == CODON_LENGTH

            and all(
                base in DNA_BASES
                for base in codon
            )
        ):

            codons.append(
                codon
            )

    return codons


def prepare_records(records):

    prepared = []

    for record in records:

        codons = (
            extract_codons(
                record[
                    "sequence"
                ]
            )
        )

        if not codons:

            raise ValueError(
                f"No valid codons extracted from "
                f"{record['name']} "
                f"({record['id']})."
            )

        prepared.append(
            {
                **record,

                "codons":
                    codons,
            }
        )

    return prepared


# ============================================================
# CONTROL CDS FEATURE EXTRACTION
# ============================================================

def clean_sequence(sequence):

    return (
        re.sub(
            r"[^A-Za-z]",
            "",
            sequence
        )
        .upper()
    )


def feature_record(
    name,
    header,
    sequence
):

    seq = (
        clean_sequence(
            sequence
        )
    )

    if (
        len(seq) < 3
        or len(seq) % 3 != 0
    ):

        return None

    codons = [
        seq[
            position:
            position + 3
        ]

        for position in range(
            0,
            len(seq),
            3
        )
    ]

    if any(
        (
            len(codon) != 3

            or any(
                base not in DNA_BASES

                for base
                in codon
            )
        )

        for codon
        in codons
    ):

        return None

    if (
        codons
        and codons[-1]
        in STOP_CODONS
    ):

        sense_codons = (
            codons[:-1]
        )

    else:

        sense_codons = (
            codons
        )

    if any(
        codon in STOP_CODONS

        for codon
        in sense_codons
    ):

        return None

    if not sense_codons:

        return None

    amino_acids = []

    for codon in sense_codons:

        amino_acid = (
            GENETIC_CODE.get(
                codon
            )
        )

        if (
            amino_acid is None
            or amino_acid == "*"
        ):

            return None

        amino_acids.append(
            amino_acid
        )

    aa_counts = (
        Counter(
            amino_acids
        )
    )

    aa_total = float(
        len(
            amino_acids
        )
    )

    aa_frequency = np.array(
        [
            aa_counts[
                amino_acid
            ]
            / aa_total

            for amino_acid
            in AA_ORDER
        ],
        dtype=float
    )

    coding_sequence = (
        "".join(
            sense_codons
        )
    )

    gc = (
        coding_sequence.count(
            "G"
        )
        + coding_sequence.count(
            "C"
        )
    ) / len(
        coding_sequence
    )

    gc3 = (
        sum(
            1

            for codon
            in sense_codons

            if codon[2]
            in "GC"
        )
        / len(
            sense_codons
        )
    )

    record_id = (
        header.split()[0]
        if header
        else ""
    )

    return {
        "name":
            name,

        "header":
            header,

        "id":
            record_id,

        "seq":
            seq,

        "codons":
            sense_codons,

        "codon_count":
            len(
                sense_codons
            ),

        "gc":
            float(
                gc
            ),

        "gc3":
            float(
                gc3
            ),

        "aa_freq":
            aa_frequency,
    }


# ============================================================
# FOCUS RECORDS FOR CONTROL SELECTION
# ============================================================

def build_focus_selection_records(
    focus_sequences
):

    records = []

    for focus in focus_sequences:

        features = (
            feature_record(
                name=
                    focus[
                        "name"
                    ],

                header=
                    focus[
                        "header"
                    ],

                sequence=
                    focus[
                        "sequence"
                    ]
            )
        )

        if features is None:

            raise ValueError(
                f"Focus sequence failed control-selection "
                f"CDS filters:\n"
                f"{focus['name']} "
                f"({focus['id']})"
            )

        records.append(
            features
        )

    if not records:

        raise ValueError(
            "No usable focus CDS records were detected."
        )

    return records


# ============================================================
# BACKGROUND CDS LOADING
# ============================================================

def parse_background_records(
    path
):

    records = []

    fasta_records = (
        read_fasta_records(
            path
        )
    )

    for (
        record_number,
        (
            header,
            sequence
        )

    ) in enumerate(
        fasta_records,
        start=1
    ):

        record_id = (
            header
            .split()[0]
        )

        record_name = (
            record_id
            if record_id
            else (
                f"background_"
                f"{record_number}"
            )
        )

        features = (
            feature_record(
                record_name,
                header,
                sequence
            )
        )

        if features is None:
            continue

        if (
            features[
                "codon_count"
            ]
            < MIN_CONTROL_CODONS
        ):

            continue

        records.append(
            features
        )

    if not records:

        raise ValueError(
            f"No usable background CDS records found in:\n"
            f"{path}"
        )

    return records


# ============================================================
# JENSEN-SHANNON DISTANCE
# ============================================================

def js_distance(
    p,
    q
):

    p = np.asarray(
        p,
        dtype=float
    )

    q = np.asarray(
        q,
        dtype=float
    )

    if (
        np.any(
            p < 0
        )
        or np.any(
            q < 0
        )
    ):

        raise ValueError(
            "Jensen-Shannon distributions "
            "cannot contain negative values."
        )

    p_total = float(
        p.sum()
    )

    q_total = float(
        q.sum()
    )

    if (
        p_total <= 0
        or q_total <= 0
    ):

        raise ValueError(
            "Jensen-Shannon distributions "
            "must have a positive total."
        )

    p = (
        p
        / p_total
    )

    q = (
        q
        / q_total
    )

    midpoint = (
        0.5
        * (
            p + q
        )
    )

    divergence = 0.0

    p_mask = (
        p > 0
    )

    if np.any(
        p_mask
    ):

        divergence += (
            0.5
            * float(
                np.sum(
                    p[
                        p_mask
                    ]
                    * np.log2(
                        p[
                            p_mask
                        ]
                        / midpoint[
                            p_mask
                        ]
                    )
                )
            )
        )

    q_mask = (
        q > 0
    )

    if np.any(
        q_mask
    ):

        divergence += (
            0.5
            * float(
                np.sum(
                    q[
                        q_mask
                    ]
                    * np.log2(
                        q[
                            q_mask
                        ]
                        / midpoint[
                            q_mask
                        ]
                    )
                )
            )
        )

    divergence = max(
        divergence,
        0.0
    )

    return float(
        math.sqrt(
            divergence
        )
    )


# ============================================================
# CONTROL MATCHING Z-SCORES
# ============================================================

def component_zscores(
    values
):

    array = np.asarray(
        values,
        dtype=float
    )

    standard_deviation = float(
        np.std(
            array
        )
    )

    if (
        standard_deviation == 0.0
        or not math.isfinite(
            standard_deviation
        )
    ):

        return (
            np.zeros_like(
                array
            )
        )

    mean = float(
        np.mean(
            array
        )
    )

    return (
        array - mean
    ) / standard_deviation


# ============================================================
# CONTROL CANDIDATE SCORING
# ============================================================

def score_candidates(
    focus,
    candidates,
    rng
):

    length_distances = []

    gc_distances = []

    gc3_distances = []

    amino_acid_distances = []

    for control in candidates:

        length_distances.append(
            abs(
                math.log(
                    control[
                        "codon_count"
                    ]
                    / focus[
                        "codon_count"
                    ]
                )
            )
        )

        gc_distances.append(
            abs(
                control[
                    "gc"
                ]
                - focus[
                    "gc"
                ]
            )
        )

        gc3_distances.append(
            abs(
                control[
                    "gc3"
                ]
                - focus[
                    "gc3"
                ]
            )
        )

        amino_acid_distances.append(
            js_distance(
                control[
                    "aa_freq"
                ],
                focus[
                    "aa_freq"
                ]
            )
        )

    length_z = (
        component_zscores(
            length_distances
        )
    )

    gc_z = (
        component_zscores(
            gc_distances
        )
    )

    gc3_z = (
        component_zscores(
            gc3_distances
        )
    )

    amino_acid_z = (
        component_zscores(
            amino_acid_distances
        )
    )

    scored = []

    for (
        index,
        control

    ) in enumerate(
        candidates
    ):

        jitter = float(
            rng.uniform(
                0.0,
                1e-9
            )
        )

        match_distance = float(
            length_z[
                index
            ]
            + gc_z[
                index
            ]
            + gc3_z[
                index
            ]
            + amino_acid_z[
                index
            ]
            + jitter
        )

        scored.append(
            (
                match_distance,
                control,
                float(
                    amino_acid_distances[
                        index
                    ]
                ),
            )
        )

    scored.sort(
        key=lambda item:
            item[0]
    )

    return scored


# ============================================================
# CONTROL SELECTION
# ============================================================

def select_controls(
    focus_records,
    background_records,
    number_of_controls,
    seed
):

    if number_of_controls < 1:

        raise ValueError(
            "Number of controls must be at least 1."
        )

    rng = np.random.default_rng(
        seed
    )

    focus_ids = {
        record[
            "id"
        ]

        for record
        in focus_records

        if record[
            "id"
        ]
    }

    focus_sequences = {
        record[
            "seq"
        ]

        for record
        in focus_records
    }

    eligible_background = [
        record

        for record
        in background_records

        if (
            record[
                "id"
            ]
            not in focus_ids

            and record[
                "seq"
            ]
            not in focus_sequences
        )
    ]

    if not eligible_background:

        raise ValueError(
            "All background records were excluded."
        )

    if (
        len(
            eligible_background
        )
        < number_of_controls
    ):

        raise ValueError(
            "The eligible background contains fewer records "
            "than the requested number of controls per focus gene."
        )

    used_ids = set()

    report_rows = []

    selected_records = []

    for (
        focus_number,
        focus

    ) in enumerate(
        focus_records,
        start=1
    ):

        print(
            f"[control selection] "
            f"{focus_number}/"
            f"{len(focus_records)}: "
            f"{focus['name']}"
        )

        candidates = [
            record

            for record
            in eligible_background

            if record[
                "id"
            ]
            not in used_ids
        ]

        # Preserve previous fallback behavior if there are not
        # enough unused controls remaining.

        if (
            len(
                candidates
            )
            < number_of_controls
        ):

            candidates = (
                eligible_background
            )

        scored = (
            score_candidates(
                focus,
                candidates,
                rng
            )
        )

        selected = (
            scored[
                :number_of_controls
            ]
        )

        for (
            control_rank,
            (
                match_distance,
                control,
                aa_js_distance
            )

        ) in enumerate(
            selected,
            start=1
        ):

            used_ids.add(
                control[
                    "id"
                ]
            )

            selected_records.append(
                {
                    "focus":
                        focus,

                    "control":
                        control,

                    "control_rank":
                        control_rank,

                    "match_distance":
                        match_distance,
                }
            )

            report_rows.append(
                {
                    "focus_name":
                        focus[
                            "name"
                        ],

                    "focus_id":
                        focus[
                            "id"
                        ],

                    "control_rank":
                        control_rank,

                    "control_id":
                        control[
                            "id"
                        ],

                    "control_header":
                        control[
                            "header"
                        ],

                    "focus_codon_count":
                        focus[
                            "codon_count"
                        ],

                    "control_codon_count":
                        control[
                            "codon_count"
                        ],

                    "abs_diff_codon_count":
                        abs(
                            control[
                                "codon_count"
                            ]
                            - focus[
                                "codon_count"
                            ]
                        ),

                    "focus_gc":
                        focus[
                            "gc"
                        ],

                    "control_gc":
                        control[
                            "gc"
                        ],

                    "abs_diff_gc":
                        abs(
                            control[
                                "gc"
                            ]
                            - focus[
                                "gc"
                            ]
                        ),

                    "focus_gc3":
                        focus[
                            "gc3"
                        ],

                    "control_gc3":
                        control[
                            "gc3"
                        ],

                    "abs_diff_gc3":
                        abs(
                            control[
                                "gc3"
                            ]
                            - focus[
                                "gc3"
                            ]
                        ),

                    "aa_js_distance":
                        aa_js_distance,

                    "match_distance":
                        match_distance,
                }
            )

    results = (
        pd.DataFrame(
            report_rows
        )
    )

    output_columns = [
        "focus_name",
        "focus_id",
        "control_rank",
        "control_id",
        "control_header",
        "focus_codon_count",
        "control_codon_count",
        "abs_diff_codon_count",
        "focus_gc",
        "control_gc",
        "abs_diff_gc",
        "focus_gc3",
        "control_gc3",
        "abs_diff_gc3",
        "aa_js_distance",
        "match_distance",
    ]

    return (
        results[
            output_columns
        ].copy(),
        selected_records
    )


# ============================================================
# PREPARE SELECTED CONTROLS FOR NULL TEST
# ============================================================

def prepare_selected_controls(
    selected_records
):

    records = []

    for item in selected_records:

        focus = (
            item[
                "focus"
            ]
        )

        control = (
            item[
                "control"
            ]
        )

        records.append(
            {
                "id":
                    control[
                        "id"
                    ],

                "name":
                    control[
                        "id"
                    ],

                "header":
                    control[
                        "header"
                    ],

                "sequence":
                    control[
                        "seq"
                    ],

                "source_file":
                    None,

                "role":
                    "control",

                "focus_id":
                    focus[
                        "id"
                    ],

                "focus_name":
                    focus[
                        "name"
                    ],

                "control_rank":
                    item[
                        "control_rank"
                    ],
            }
        )

    return (
        prepare_records(
            records
        )
    )


# ============================================================
# tRNA INPUT
# ============================================================

TRNA_PAIR_PATTERN = re.compile(
    r"""
    (?<![A-Za-z])
    ([ACGTUacgtu]{3})
    (?![A-Za-z])
    \s*
    [,;:\t|]?
    \s*
    ([+-]?\d+(?:\.0+)?)
    """,
    re.VERBOSE
)


def read_trna_table(path):

    text = (
        path.read_text(
            encoding="utf-8-sig"
        )
    )

    matches = list(
        TRNA_PAIR_PATTERN.finditer(
            text
        )
    )

    if not matches:

        raise ValueError(
            "No codon/count pairs detected "
            "in the tRNA file."
        )

    counts = {
        codon: 0

        for codon
        in ALL_CODONS
    }

    observed = set()

    duplicates = set()

    for match in matches:

        codon = (
            match.group(1)
            .upper()
            .replace(
                "U",
                "T"
            )
        )

        numeric_count = float(
            match.group(2)
        )

        if not (
            numeric_count.is_integer()
        ):

            raise ValueError(
                f"Non-integer tRNA copy number "
                f"for {codon}: "
                f"{numeric_count}"
            )

        copies = int(
            numeric_count
        )

        if copies < 0:

            raise ValueError(
                f"Negative tRNA copy number "
                f"for {codon}: "
                f"{copies}"
            )

        if codon in observed:

            duplicates.add(
                codon
            )

        counts[
            codon
        ] += copies

        observed.add(
            codon
        )

    return {
        "counts":
            counts,

        "pairs_read":
            len(
                matches
            ),

        "unique_codons":
            len(
                observed
            ),

        "total_copies":
            sum(
                counts.values()
            ),

        "duplicates":
            sorted(
                duplicates
            ),
    }


# ============================================================
# sij INPUT
# ============================================================

SIJ_PATTERN = re.compile(
    r"""
    ^\s*
    [*-]?
    \s*
    s[\s_-]*
    (GU|UG|IC|IA)
    \s*
    (?:
        :
        |
        =
        |
        \s+
    )
    \s*
    (
        [+-]?
        (?:
            \d+(?:\.\d*)?
            |
            \.\d+
        )
        (?:[eE][+-]?\d+)?
    )
    \s*$
    """,
    re.VERBOSE
    | re.IGNORECASE
)


def read_sij_table(path):

    text = (
        path.read_text(
            encoding="utf-8-sig"
        )
    )

    final_marker = re.search(
        r"""
        final
        \s+
        refined
        \s+
        combination
        \s*:
        """,
        text,
        flags=(
            re.IGNORECASE
            | re.VERBOSE
        )
    )

    if final_marker:

        text = text[
            final_marker.end():
        ]

    penalties = {}

    for (
        line_number,
        line

    ) in enumerate(
        text.splitlines(),
        start=1
    ):

        match = (
            SIJ_PATTERN.match(
                line
            )
        )

        if match is None:
            continue

        parameter = (
            f"s_"
            f"{match.group(1).upper()}"
        )

        value = float(
            match.group(2)
        )

        if not math.isfinite(
            value
        ):

            raise ValueError(
                f"Non-finite value for "
                f"{parameter}."
            )

        if not (
            0.0
            <= value
            <= 0.9999
        ):

            raise ValueError(
                f"{parameter} must be between "
                f"0 and 0.9999. "
                f"Found: {value}"
            )

        if parameter in penalties:

            raise ValueError(
                f"Duplicate value for "
                f"{parameter} near line "
                f"{line_number}."
            )

        penalties[
            parameter
        ] = value

    missing = [
        parameter

        for parameter
        in SIJ_PARAMETERS

        if parameter
        not in penalties
    ]

    if missing:

        raise ValueError(
            "sij file is missing: "
            + ", ".join(
                missing
            )
        )

    return penalties


# ============================================================
# tRNA AVAILABILITY MODEL
# ============================================================

def anticodon_rna_for_codon(codon):

    complement = {
        "A": "U",
        "T": "A",
        "C": "G",
        "G": "C",
    }

    return "".join(
        complement[
            base
        ]

        for base
        in codon
    )[::-1]


def build_wobble_edges(
    trna_counts
):

    edges = []

    for source_codon in SENSE_CODONS:

        copies = float(
            trna_counts.get(
                source_codon,
                0
            )
        )

        if copies <= 0:
            continue

        anticodon = (
            anticodon_rna_for_codon(
                source_codon
            )
        )

        amino_acid = (
            GENETIC_CODE[
                source_codon
            ]
        )

        source_rna = (
            source_codon.replace(
                "T",
                "U"
            )
        )

        if (
            anticodon
            in I34_CANDIDATE_ANTICODONS_RNA
        ):

            effective_bases = (
                "I",
            )

        else:

            effective_bases = (
                anticodon[0],
            )

        for base_34 in effective_bases:

            for codon_third_base in (
                WOBBLE_THIRD_BASES[
                    base_34
                ]
            ):

                target_rna = (
                    source_rna[:2]
                    + codon_third_base
                )

                target_codon = (
                    target_rna.replace(
                        "U",
                        "T"
                    )
                )

                if (
                    target_codon
                    not in SENSE_CODONS
                ):

                    continue

                if (
                    GENETIC_CODE[
                        target_codon
                    ]
                    != amino_acid
                ):

                    continue

                pair_name = (
                    f"{base_34}:"
                    f"{codon_third_base}"
                )

                edges.append(
                    {
                        "codon":
                            target_codon,

                        "copies":
                            copies,

                        "parameter":
                            PAIR_TO_SIJ_PARAMETER.get(
                                pair_name
                            ),
                    }
                )

    return edges


def calculate_relative_codon_weights(
    trna_counts,
    sij_values
):

    edges = (
        build_wobble_edges(
            trna_counts
        )
    )

    raw_weights = {
        codon: 0.0

        for codon
        in SENSE_CODONS
    }

    for edge in edges:

        parameter = (
            edge[
                "parameter"
            ]
        )

        if parameter is None:

            sij = 0.0

        else:

            sij = (
                sij_values[
                    parameter
                ]
            )

        raw_weights[
            edge[
                "codon"
            ]
        ] += (
            edge[
                "copies"
            ]
            * (
                1.0
                - sij
            )
        )

    zero_codons = [
        codon

        for codon, value
        in raw_weights.items()

        if value <= 0
    ]

    if zero_codons:

        raise ValueError(
            "Zero tRNA availability detected for: "
            + ", ".join(
                zero_codons
            )
        )

    maximum = max(
        raw_weights.values()
    )

    relative_weights = {
        codon:
            raw_weights[
                codon
            ]
            / maximum

        for codon
        in SENSE_CODONS
    }

    inverse_weights = {
        codon:
            1.0
            / relative_weights[
                codon
            ]

        for codon
        in SENSE_CODONS
    }

    return (
        raw_weights,
        relative_weights,
        inverse_weights,
        edges
    )


# ============================================================
# SYNONYMOUS CODON FAMILIES
# ============================================================

def build_codon_families():

    families = {}

    for codon, amino_acid in (
        GENETIC_CODE.items()
    ):

        if amino_acid == "*":
            continue

        families.setdefault(
            amino_acid,
            []
        ).append(
            codon
        )

    for amino_acid in families:

        families[
            amino_acid
        ] = sorted(
            families[
                amino_acid
            ]
        )

    return families


CODON_FAMILIES = (
    build_codon_families()
)


# ============================================================
# HOST-GLOBAL SYNONYMOUS DISTRIBUTION
# ============================================================

def build_host_global_probabilities(
    hosts
):

    counts_by_amino_acid = {
        amino_acid:
            Counter()

        for amino_acid
        in CODON_FAMILIES
    }

    for host in hosts:

        for codon in (
            host[
                "codons"
            ]
        ):

            amino_acid = (
                GENETIC_CODE.get(
                    codon
                )
            )

            if (
                amino_acid
                in counts_by_amino_acid
            ):

                counts_by_amino_acid[
                    amino_acid
                ][
                    codon
                ] += 1

    probabilities = {}

    for (
        amino_acid,
        codons

    ) in CODON_FAMILIES.items():

        raw = np.array(
            [
                (
                    counts_by_amino_acid[
                        amino_acid
                    ][
                        codon
                    ]
                    + HOST_GLOBAL_PSEUDOCOUNT
                )

                for codon
                in codons
            ],
            dtype=float
        )

        total = float(
            raw.sum()
        )

        if total <= 0:

            raise ValueError(
                f"Unable to construct synonymous "
                f"probabilities for amino acid "
                f"{amino_acid}."
            )

        probabilities[
            amino_acid
        ] = (
            codons,
            raw
            / total
        )

    return probabilities


# ============================================================
# SYNONYMOUS RECODING
# ============================================================

def synonymous_recode(
    codons,
    probabilities,
    rng
):

    recoded = []

    for codon in codons:

        amino_acid = (
            GENETIC_CODE.get(
                codon
            )
        )

        if (
            amino_acid is None
            or amino_acid == "*"
        ):

            recoded.append(
                codon
            )

            continue

        family, probability = (
            probabilities[
                amino_acid
            ]
        )

        if len(
            family
        ) <= 1:

            recoded.append(
                codon
            )

        else:

            recoded.append(
                str(
                    rng.choice(
                        family,
                        p=probability
                    )
                )
            )

    return recoded


# ============================================================
# FCGR CONSTRUCTION
# ============================================================

def create_relative_frequency_fcgr(
    codons,
    inverse_weights=None
):

    matrix = np.zeros(
        (
            MATRIX_SIZE,
            MATRIX_SIZE
        ),
        dtype=float
    )

    for codon in codons:

        coordinate = (
            CODON_TO_COORDINATE.get(
                codon
            )
        )

        if coordinate is None:
            continue

        x, y = coordinate

        if inverse_weights is None:

            contribution = 1.0

        elif codon in SENSE_CODONS:

            contribution = (
                inverse_weights[
                    codon
                ]
            )

        else:

            contribution = 1.0

        matrix[
            y,
            x
        ] += contribution

    total = float(
        matrix.sum()
    )

    if total > 0:

        matrix /= total

    return matrix


# ============================================================
# DIRECT COMPARISON METRICS
# ============================================================

def calculate_direct_metrics(
    virus_matrix,
    host_matrix
):
    """
    Calculate the standardized Virus-Host comparison metrics.

    Spearman:
        scipy.stats.spearmanr

    Pearson:
        scipy.stats.pearsonr

    SSIM:
        skimage.metrics.structural_similarity

    SSIM requires finite input matrices. Non-finite values are
    treated as invalid input and stop execution rather than being
    ignored or silently converted to NaN.

    The SSIM data range is defined from the Host matrix:

        max(Host) - min(Host)

    If and only if this finite range is exactly zero, a value of
    1.0 is used because SSIM requires a non-zero data_range.
    """

    virus_matrix = np.asarray(
        virus_matrix,
        dtype=float
    )

    host_matrix = np.asarray(
        host_matrix,
        dtype=float
    )

    virus_flat = (
        virus_matrix.flatten()
    )

    host_flat = (
        host_matrix.flatten()
    )

    # --------------------------------------------------------
    # Spearman
    # --------------------------------------------------------

    with warnings.catch_warnings():

        warnings.simplefilter(
            "ignore"
        )

        spearman_result = (
            stats.spearmanr(
                virus_flat,
                host_flat
            )
        )

    spearman_r = float(
        spearman_result.statistic
    )

    # --------------------------------------------------------
    # Pearson
    # --------------------------------------------------------

    with warnings.catch_warnings():

        warnings.simplefilter(
            "ignore"
        )

        pearson_result = (
            stats.pearsonr(
                virus_flat,
                host_flat
            )
        )

    pearson_r = float(
        pearson_result.statistic
    )

    # --------------------------------------------------------
    # SSIM input validation
    # --------------------------------------------------------

    if not np.all(
        np.isfinite(
            virus_matrix
        )
    ):

        raise ValueError(
            "Virus matrix contains non-finite values; "
            "SSIM cannot be calculated."
        )

    if not np.all(
        np.isfinite(
            host_matrix
        )
    ):

        raise ValueError(
            "Host matrix contains non-finite values; "
            "SSIM cannot be calculated."
        )

    data_range = float(
        np.max(
            host_matrix
        )
        - np.min(
            host_matrix
        )
    )

    if data_range == 0.0:

        data_range = 1.0

    ssim_value = float(
        ssim(
            virus_matrix,
            host_matrix,
            data_range=data_range
        )
    )

    return {
        "spearman_r":
            spearman_r,

        "pearson_r":
            pearson_r,

        "ssim":
            ssim_value,
    }


# ============================================================
# COMPETITION LAYOUT
# ============================================================

def build_competition_layout(
    relative_weights
):

    weight_matrix = np.ones(
        (
            MATRIX_SIZE,
            MATRIX_SIZE
        ),
        dtype=float
    )

    sense_mask = np.zeros(
        (
            MATRIX_SIZE,
            MATRIX_SIZE
        ),
        dtype=bool
    )

    for codon in SENSE_CODONS:

        x, y = (
            CODON_TO_COORDINATE[
                codon
            ]
        )

        weight_matrix[
            y,
            x
        ] = (
            relative_weights[
                codon
            ]
        )

        sense_mask[
            y,
            x
        ] = True

    return (
        weight_matrix,
        sense_mask
    )


# ============================================================
# SENSE-CODON USAGE
# ============================================================

def normalize_sense_usage(
    matrix,
    sense_mask
):

    values = np.asarray(
        matrix,
        dtype=float
    ).copy()

    values[
        ~sense_mask
    ] = 0.0

    total = float(
        values.sum()
    )

    if total > 0:

        values /= total

    return values


# ============================================================
# COMPETITION SCORE
# ============================================================

def competition_normalized_score(
    virus_unweighted_matrix,
    host_unweighted_matrix,
    relative_weights,
    weight_matrix,
    sense_mask
):

    virus_frequency = (
        normalize_sense_usage(
            virus_unweighted_matrix,
            sense_mask
        )
    )

    host_frequency = (
        normalize_sense_usage(
            host_unweighted_matrix,
            sense_mask
        )
    )

    minimum_weight = min(
        relative_weights.values()
    )

    theoretical_maximum = (
        1.0
        / minimum_weight
    )

    scarcity = (
        1.0
        / weight_matrix
    )

    raw = (
        host_frequency
        * virus_frequency
        * scarcity
    )

    normalized = (
        raw
        / theoretical_maximum
    )

    return float(
        normalized.sum()
    )


# ============================================================
# VIRUS ARTIFACT
# ============================================================

def build_virus_artifact(
    virus,
    mode,
    inverse_weights
):

    unweighted_matrix = (
        create_relative_frequency_fcgr(
            virus[
                "codons"
            ]
        )
    )

    if mode == "penalized":

        comparison_matrix = (
            create_relative_frequency_fcgr(
                virus[
                    "codons"
                ],
                inverse_weights=
                    inverse_weights
            )
        )

    else:

        comparison_matrix = (
            unweighted_matrix
        )

    return {
        "comparison_matrix":
            comparison_matrix,

        "unweighted_matrix":
            unweighted_matrix,
    }


# ============================================================
# HOST ARTIFACT
# ============================================================

def build_host_artifact(
    host_codons
):

    matrix = (
        create_relative_frequency_fcgr(
            host_codons
        )
    )

    return {
        "comparison_matrix":
            matrix,

        "unweighted_matrix":
            matrix,
    }


# ============================================================
# COMPLETE METRICS
# ============================================================

def compute_metrics(
    virus_artifact,
    host_artifact,
    relative_weights,
    weight_matrix,
    sense_mask
):

    metrics = (
        calculate_direct_metrics(
            virus_artifact[
                "comparison_matrix"
            ],
            host_artifact[
                "comparison_matrix"
            ]
        )
    )

    metrics[
        "competition_normalized_score"
    ] = (
        competition_normalized_score(
            virus_artifact[
                "unweighted_matrix"
            ],
            host_artifact[
                "unweighted_matrix"
            ],
            relative_weights,
            weight_matrix,
            sense_mask
        )
    )

    return metrics


# ============================================================
# EMPIRICAL P-VALUE
# ============================================================

def empirical_p(
    observed,
    null_values
):

    null = np.asarray(
        null_values,
        dtype=float
    )

    null = (
        null[
            np.isfinite(
                null
            )
        ]
    )

    if (
        not math.isfinite(
            float(
                observed
            )
        )
        or null.size == 0
    ):

        return (
            float(
                "nan"
            ),
            0
        )

    b_extreme = int(
        np.sum(
            null
            >= observed
        )
    )

    p_perm = (
        b_extreme
        + 1.0
    ) / (
        null.size
        + 1.0
    )

    return (
        float(
            p_perm
        ),
        b_extreme
    )


# ============================================================
# NULL DISTRIBUTION SUMMARY
# ============================================================

def null_summary(values):

    null = np.asarray(
        values,
        dtype=float
    )

    null = (
        null[
            np.isfinite(
                null
            )
        ]
    )

    if null.size == 0:

        return {
            "null_mean":
                np.nan,

            "null_sd":
                np.nan,

            "null_q025":
                np.nan,

            "null_median":
                np.nan,

            "null_q975":
                np.nan,
        }

    return {
        "null_mean":
            float(
                np.mean(
                    null
                )
            ),

        "null_sd":
            (
                float(
                    np.std(
                        null,
                        ddof=1
                    )
                )

                if null.size > 1

                else 0.0
            ),

        "null_q025":
            float(
                np.quantile(
                    null,
                    0.025
                )
            ),

        "null_median":
            float(
                np.median(
                    null
                )
            ),

        "null_q975":
            float(
                np.quantile(
                    null,
                    0.975
                )
            ),
    }


# ============================================================
# BENJAMINI-HOCHBERG
# ============================================================

def benjamini_hochberg(
    values
):

    values = np.asarray(
        values,
        dtype=float
    )

    corrected = np.full(
        len(
            values
        ),
        np.nan,
        dtype=float
    )

    valid = (
        np.isfinite(
            values
        )
        & (
            values >= 0
        )
        & (
            values <= 1
        )
    )

    if np.any(
        valid
    ):

        corrected[
            valid
        ] = (
            multipletests(
                values[
                    valid
                ],
                method="fdr_bh"
            )[1]
        )

    return corrected


def add_q_values(
    results
):

    results = (
        results.copy()
    )

    results[
        "q_perm"
    ] = np.nan

    groups = (
        results.groupby(
            [
                "mode",
                "metric",
            ],
            sort=False
        ).groups
    )

    for (
        _group,
        indexes

    ) in groups.items():

        results.loc[
            indexes,
            "q_perm"
        ] = (
            benjamini_hochberg(
                pd.to_numeric(
                    results.loc[
                        indexes,
                        "p_perm"
                    ],
                    errors="coerce"
                ).to_numpy(
                    dtype=float
                )
            )
        )

    return results


# ============================================================
# NULL TEST
# ============================================================

def run_null_test(
    viruses,
    hosts,
    relative_weights,
    inverse_weights,
    synonymous_probabilities,
    seed,
    permutations,
    dataset_label
):

    rng = (
        np.random.default_rng(
            seed
        )
    )

    (
        competition_weight_matrix,
        sense_mask

    ) = (
        build_competition_layout(
            relative_weights
        )
    )

    result_rows = []

    total_comparisons = (
        len(
            viruses
        )
        * len(
            hosts
        )
    )

    for mode in MODES:

        print()

        print(
            "=" * 78
        )

        print(
            f"{dataset_label.upper()} NULL TEST: "
            f"{mode.upper()}"
        )

        print(
            "=" * 78
        )

        print(
            f"Sequences                       : "
            f"{len(hosts)}"
        )

        print(
            f"Comparisons                     : "
            f"{total_comparisons}"
        )

        print(
            f"Permutations per sequence       : "
            f"{permutations}"
        )

        print()

        virus_artifacts = {
            virus[
                "name"
            ]:
                build_virus_artifact(
                    virus,
                    mode,
                    inverse_weights
                )

            for virus
            in viruses
        }

        for (
            host_number,
            host

        ) in enumerate(
            hosts,
            start=1
        ):

            print(
                f"[{dataset_label} | {mode}] "
                f"{host_number}/"
                f"{len(hosts)}: "
                f"{host['name']}"
            )

            observed_host_artifact = (
                build_host_artifact(
                    host[
                        "codons"
                    ]
                )
            )

            observed_by_virus = {}

            for virus in viruses:

                virus_name = (
                    virus[
                        "name"
                    ]
                )

                observed_by_virus[
                    virus_name
                ] = (
                    compute_metrics(
                        virus_artifacts[
                            virus_name
                        ],
                        observed_host_artifact,
                        relative_weights,
                        competition_weight_matrix,
                        sense_mask
                    )
                )

            null_values = {
                (
                    virus[
                        "name"
                    ],
                    metric
                ):
                    []

                for virus
                in viruses

                for metric
                in METRICS
            }

            for _permutation_number in range(
                1,
                permutations + 1
            ):

                recoded_codons = (
                    synonymous_recode(
                        host[
                            "codons"
                        ],
                        synonymous_probabilities,
                        rng
                    )
                )

                null_host_artifact = (
                    build_host_artifact(
                        recoded_codons
                    )
                )

                for virus in viruses:

                    virus_name = (
                        virus[
                            "name"
                        ]
                    )

                    permutation_metrics = (
                        compute_metrics(
                            virus_artifacts[
                                virus_name
                            ],
                            null_host_artifact,
                            relative_weights,
                            competition_weight_matrix,
                            sense_mask
                        )
                    )

                    for metric in METRICS:

                        null_values[
                            (
                                virus_name,
                                metric
                            )
                        ].append(
                            permutation_metrics[
                                metric
                            ]
                        )

            for virus in viruses:

                virus_name = (
                    virus[
                        "name"
                    ]
                )

                for metric in METRICS:

                    observed = (
                        observed_by_virus[
                            virus_name
                        ][
                            metric
                        ]
                    )

                    values = (
                        null_values[
                            (
                                virus_name,
                                metric
                            )
                        ]
                    )

                    (
                        p_perm,
                        b_extreme

                    ) = (
                        empirical_p(
                            observed,
                            values
                        )
                    )

                    summary = (
                        null_summary(
                            values
                        )
                    )

                    row = {
                        "mode":
                            mode,

                        "virus_id":
                            virus[
                                "id"
                            ],

                        "virus_name":
                            virus_name,

                        "metric":
                            metric,

                        "observed":
                            observed,

                        "null_mean":
                            summary[
                                "null_mean"
                            ],

                        "null_sd":
                            summary[
                                "null_sd"
                            ],

                        "null_q025":
                            summary[
                                "null_q025"
                            ],

                        "null_median":
                            summary[
                                "null_median"
                            ],

                        "null_q975":
                            summary[
                                "null_q975"
                            ],

                        "b_extreme":
                            b_extreme,

                        "p_perm":
                            p_perm,
                    }

                    if dataset_label == "focus":

                        row[
                            "focus_id"
                        ] = (
                            host[
                                "id"
                            ]
                        )

                        row[
                            "focus_name"
                        ] = (
                            host[
                                "name"
                            ]
                        )

                    elif dataset_label == "control":

                        row[
                            "focus_id"
                        ] = (
                            host[
                                "focus_id"
                            ]
                        )

                        row[
                            "focus_name"
                        ] = (
                            host[
                                "focus_name"
                            ]
                        )

                        row[
                            "control_id"
                        ] = (
                            host[
                                "id"
                            ]
                        )

                        row[
                            "control_rank"
                        ] = (
                            host[
                                "control_rank"
                            ]
                        )

                    else:

                        raise ValueError(
                            f"Unknown dataset label: "
                            f"{dataset_label}"
                        )

                    result_rows.append(
                        row
                    )

    results = (
        pd.DataFrame(
            result_rows
        )
    )

    results = (
        add_q_values(
            results
        )
    )

    if dataset_label == "focus":

        final_columns = [
            "mode",
            "virus_id",
            "virus_name",
            "focus_id",
            "focus_name",
            "metric",
            "observed",
            "null_mean",
            "null_sd",
            "null_q025",
            "null_median",
            "null_q975",
            "b_extreme",
            "p_perm",
            "q_perm",
        ]

    else:

        final_columns = [
            "mode",
            "virus_id",
            "virus_name",
            "focus_id",
            "focus_name",
            "control_id",
            "control_rank",
            "metric",
            "observed",
            "null_mean",
            "null_sd",
            "null_q025",
            "null_median",
            "null_q975",
            "b_extreme",
            "p_perm",
            "q_perm",
        ]

    return (
        results[
            final_columns
        ].copy()
    )


# ============================================================
# BUILD CONTROL COMPARISON TABLE
# ============================================================

def build_control_comparison(
    focus_results,
    control_results
):

    focus = (
        focus_results.copy()
    )

    controls = (
        control_results.copy()
    )

    focus = focus.rename(
        columns={
            "observed":
                "focus_observed",

            "null_mean":
                "focus_null_mean",

            "null_median":
                "focus_null_median",
        }
    )

    focus_columns = [
        "mode",
        "virus_id",
        "virus_name",
        "focus_id",
        "focus_name",
        "metric",
        "focus_observed",
        "focus_null_mean",
        "focus_null_median",
    ]

    focus = (
        focus[
            focus_columns
        ].copy()
    )

    controls = controls.rename(
        columns={
            "observed":
                "control_observed",

            "null_mean":
                "control_null_mean",

            "null_sd":
                "control_null_sd",

            "null_q025":
                "control_null_q025",

            "null_median":
                "control_null_median",

            "null_q975":
                "control_null_q975",

            "b_extreme":
                "control_b_extreme",

            "p_perm":
                "control_p_perm",

            "q_perm":
                "control_q_perm",
        }
    )

    comparison = controls.merge(
        focus,
        on=[
            "mode",
            "virus_id",
            "virus_name",
            "focus_id",
            "focus_name",
            "metric",
        ],
        how="left",
        validate="many_to_one"
    )

    missing_focus = int(
        comparison[
            "focus_observed"
        ].isna().sum()
    )

    if missing_focus > 0:

        raise ValueError(
            f"{missing_focus} control rows could not "
            f"be matched to focus null-test results."
        )

    comparison[
        "abs_observed_difference"
    ] = np.abs(
        comparison[
            "focus_observed"
        ]
        - comparison[
            "control_observed"
        ]
    )

    comparison[
        "abs_null_mean_difference"
    ] = np.abs(
        comparison[
            "focus_null_mean"
        ]
        - comparison[
            "control_null_mean"
        ]
    )

    comparison[
        "abs_null_median_difference"
    ] = np.abs(
        comparison[
            "focus_null_median"
        ]
        - comparison[
            "control_null_median"
        ]
    )

    final_columns = [
        "mode",
        "virus_id",
        "virus_name",
        "focus_id",
        "focus_name",
        "control_id",
        "control_rank",
        "metric",

        "control_observed",
        "focus_observed",
        "abs_observed_difference",

        "control_null_mean",
        "focus_null_mean",
        "abs_null_mean_difference",

        "control_null_median",
        "focus_null_median",
        "abs_null_median_difference",

        "control_null_sd",
        "control_null_q025",
        "control_null_q975",
        "control_b_extreme",
        "control_p_perm",
        "control_q_perm",
    ]

    return (
        comparison[
            final_columns
        ].copy()
    )


# ============================================================
# ARGUMENTS
# ============================================================

def parse_arguments():

    # -h is used for Host.
    # Help is therefore --help.

    parser = argparse.ArgumentParser(
        description=(
            "Integrated synonymous-recoding null test, "
            "matched-control selection, and individual "
            "focus-control comparison."
        ),
        add_help=False
    )

    parser.add_argument(
        "--help",
        action="help",
        help=(
            "Show this help message and exit."
        )
    )

    parser.add_argument(
        "-v",
        "--virus",
        required=True,
        help=(
            "Virus FASTA file or directory."
        )
    )

    parser.add_argument(
        "-h",
        "--host",
        required=True,
        help=(
            "Focus Host FASTA file or directory."
        )
    )

    parser.add_argument(
        "-t",
        "--trna",
        required=True,
        help=(
            "Plain-text tRNA codon/copy-number file."
        )
    )

    parser.add_argument(
        "-j",
        "--sij",
        required=True,
        help=(
            "Plain-text file containing "
            "s_GU, s_UG, s_IC and s_IA."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=(
            "Background CDS multi-FASTA file "
            "for matched-control selection."
        )
    )

    parser.add_argument(
        "-s",
        "--seed",
        required=True,
        type=int,
        help=(
            "Random seed."
        )
    )

    parser.add_argument(
        "-p",
        "--perm",
        required=True,
        type=int,
        help=(
            "Number of synonymous-recoding "
            "permutations per sequence."
        )
    )

    parser.add_argument(
        "-c",
        "--controls",
        required=True,
        type=int,
        help=(
            "Number of matched controls "
            "selected per focus gene."
        )
    )

    args = (
        parser.parse_args()
    )

    if args.seed < 0:

        parser.error(
            "--seed must be a non-negative integer."
        )

    if args.perm < 1:

        parser.error(
            "--perm must be at least 1."
        )

    if args.controls < 1:

        parser.error(
            "--controls must be at least 1."
        )

    return args


# ============================================================
# MAIN
# ============================================================

def main():

    args = (
        parse_arguments()
    )

    start_time = (
        time.time()
    )

    # ========================================================
    # PATHS
    # ========================================================

    virus_source = (
        resolve_path(
            args.virus
        )
    )

    host_source = (
        resolve_path(
            args.host
        )
    )

    trna_path = (
        resolve_path(
            args.trna
        )
    )

    sij_path = (
        resolve_path(
            args.sij
        )
    )

    background_path = (
        resolve_path(
            args.input
        )
    )

    null_workbook_path = (
        script_directory()
        / NULL_WORKBOOK_NAME
    )

    controls_workbook_path = (
        script_directory()
        / CONTROLS_WORKBOOK_NAME
    )

    # ========================================================
    # VALIDATION
    # ========================================================

    validate_sequence_source(
        virus_source,
        "Virus"
    )

    validate_sequence_source(
        host_source,
        "Host"
    )

    validate_file(
        trna_path,
        "tRNA"
    )

    validate_file(
        sij_path,
        "sij"
    )

    validate_file(
        background_path,
        "Background CDS"
    )

    # ========================================================
    # LOAD VIRUSES AND FOCUS HOSTS
    # ========================================================

    virus_sequences = (
        load_sequence_source(
            virus_source,
            role="virus"
        )
    )

    focus_sequences = (
        load_sequence_source(
            host_source,
            role="focus"
        )
    )

    viruses = (
        prepare_records(
            virus_sequences
        )
    )

    focus_hosts = (
        prepare_records(
            focus_sequences
        )
    )

    focus_selection_records = (
        build_focus_selection_records(
            focus_sequences
        )
    )

    if (
        len(
            focus_hosts
        )
        != len(
            focus_selection_records
        )
    ):

        raise RuntimeError(
            "Focus sequence count mismatch between "
            "null-test and control-selection preparation."
        )

    background_records = (
        parse_background_records(
            background_path
        )
    )

    # ========================================================
    # LOAD tRNA AND sij
    # ========================================================

    trna_result = (
        read_trna_table(
            trna_path
        )
    )

    sij_values = (
        read_sij_table(
            sij_path
        )
    )

    (
        _raw_weights,
        relative_weights,
        inverse_weights,
        wobble_edges

    ) = (
        calculate_relative_codon_weights(
            trna_result[
                "counts"
            ],
            sij_values
        )
    )

    # ========================================================
    # FOCUS HOST-GLOBAL MODEL
    # ========================================================

    focus_synonymous_probabilities = (
        build_host_global_probabilities(
            focus_hosts
        )
    )

    # ========================================================
    # RUN SUMMARY
    # ========================================================

    print()

    print(
        "=" * 78
    )

    print(
        "INTEGRATED NULL TEST AND MATCHED-CONTROL ANALYSIS"
    )

    print(
        "=" * 78
    )

    print(
        f"Virus sequences                 : "
        f"{len(viruses)}"
    )

    print(
        f"Focus sequences                 : "
        f"{len(focus_hosts)}"
    )

    print(
        f"Usable background CDSs          : "
        f"{len(background_records)}"
    )

    print(
        f"Controls per focus gene         : "
        f"{args.controls}"
    )

    print(
        f"Expected selected controls      : "
        f"{len(focus_hosts) * args.controls}"
    )

    print(
        f"Seed                            : "
        f"{args.seed}"
    )

    print(
        f"Permutations                    : "
        f"{args.perm}"
    )

    print(
        f"Host-global pseudocount         : "
        f"{HOST_GLOBAL_PSEUDOCOUNT}"
    )

    print(
        f"tRNA copies                     : "
        f"{trna_result['total_copies']}"
    )

    print(
        f"Positive decoding edges         : "
        f"{len(wobble_edges)}"
    )

    for parameter in SIJ_PARAMETERS:

        print(
            f"{parameter:<32}: "
            f"{sij_values[parameter]:.6f}"
        )

    print(
        f"Metrics                         : "
        + ", ".join(
            METRICS
        )
    )

    print(
        f"Output null test                : "
        f"{null_workbook_path}"
    )

    print(
        f"Output controls                 : "
        f"{controls_workbook_path}"
    )

    print(
        "=" * 78
    )

    # ========================================================
    # 1. FOCUS NULL TEST
    # ========================================================

    focus_null_results = (
        run_null_test(
            viruses=
                viruses,

            hosts=
                focus_hosts,

            relative_weights=
                relative_weights,

            inverse_weights=
                inverse_weights,

            synonymous_probabilities=
                focus_synonymous_probabilities,

            seed=
                args.seed,

            permutations=
                args.perm,

            dataset_label=
                "focus"
        )
    )

    expected_focus_rows = (
        len(
            MODES
        )
        * len(
            viruses
        )
        * len(
            focus_hosts
        )
        * len(
            METRICS
        )
    )

    if (
        len(
            focus_null_results
        )
        != expected_focus_rows
    ):

        raise RuntimeError(
            f"Unexpected focus null-test row count. "
            f"Expected {expected_focus_rows}, "
            f"obtained {len(focus_null_results)}."
        )

    # ========================================================
    # 2. CONTROL SELECTION
    # ========================================================

    print()

    print(
        "=" * 78
    )

    print(
        "SELECTING MATCHED CONTROLS"
    )

    print(
        "=" * 78
    )

    (
        control_selection_results,
        selected_records

    ) = (
        select_controls(
            focus_records=
                focus_selection_records,

            background_records=
                background_records,

            number_of_controls=
                args.controls,

            seed=
                args.seed
        )
    )

    expected_controls = (
        len(
            focus_selection_records
        )
        * args.controls
    )

    if (
        len(
            control_selection_results
        )
        != expected_controls
    ):

        raise RuntimeError(
            f"Unexpected selected-control count. "
            f"Expected {expected_controls}, "
            f"obtained "
            f"{len(control_selection_results)}."
        )

    control_hosts = (
        prepare_selected_controls(
            selected_records
        )
    )

    # ========================================================
    # 3. CONTROL HOST-GLOBAL MODEL
    # ========================================================

    control_synonymous_probabilities = (
        build_host_global_probabilities(
            control_hosts
        )
    )

    # ========================================================
    # 4. CONTROL NULL TEST
    # ========================================================

    control_null_results = (
        run_null_test(
            viruses=
                viruses,

            hosts=
                control_hosts,

            relative_weights=
                relative_weights,

            inverse_weights=
                inverse_weights,

            synonymous_probabilities=
                control_synonymous_probabilities,

            seed=
                args.seed,

            permutations=
                args.perm,

            dataset_label=
                "control"
        )
    )

    expected_control_rows = (
        len(
            MODES
        )
        * len(
            viruses
        )
        * len(
            control_hosts
        )
        * len(
            METRICS
        )
    )

    if (
        len(
            control_null_results
        )
        != expected_control_rows
    ):

        raise RuntimeError(
            f"Unexpected control null-test row count. "
            f"Expected {expected_control_rows}, "
            f"obtained "
            f"{len(control_null_results)}."
        )

    # ========================================================
    # 5. BUILD INDIVIDUAL FOCUS-CONTROL COMPARISON
    # ========================================================

    comparison_results = (
        build_control_comparison(
            focus_results=
                focus_null_results,

            control_results=
                control_null_results
        )
    )

    if (
        len(
            comparison_results
        )
        != expected_control_rows
    ):

        raise RuntimeError(
            f"Unexpected comparison row count. "
            f"Expected {expected_control_rows}, "
            f"obtained "
            f"{len(comparison_results)}."
        )

    # ========================================================
    # 6. WRITE CONTROL SELECTION
    # ========================================================

    control_selection_results.to_excel(
        controls_workbook_path,
        sheet_name=
            CONTROL_SELECTION_SHEET,
        index=False,
        freeze_panes=(
            1,
            0
        )
    )

    # ========================================================
    # 7. WRITE NULL TEST RESULTS
    # ========================================================

    with pd.ExcelWriter(
        null_workbook_path
    ) as writer:

        focus_null_results.to_excel(
            writer,
            sheet_name="Focus",
            index=False,
            freeze_panes=(
                1,
                0
            )
        )

        comparison_results.to_excel(
            writer,
            sheet_name="Control",
            index=False,
            freeze_panes=(
                1,
                0
            )
        )

    # ========================================================
    # FINAL REPORT
    # ========================================================

    elapsed_seconds = (
        time.time()
        - start_time
    )

    print()

    print(
        "=" * 78
    )

    print(
        "ANALYSIS COMPLETED"
    )

    print(
        "=" * 78
    )

    print(
        f"Focus sequences                 : "
        f"{len(focus_hosts)}"
    )

    print(
        f"Selected controls               : "
        f"{len(control_hosts)}"
    )

    print(
        f"Unique control IDs              : "
        f"{control_selection_results['control_id'].nunique()}"
    )

    print(
        f"Focus result rows               : "
        f"{len(focus_null_results)}"
    )

    print(
        f"Control comparison rows         : "
        f"{len(comparison_results)}"
    )

    print(
        f"null_test.xlsx                  : "
        f"{null_workbook_path}"
    )

    print(
        f"controls.xlsx                   : "
        f"{controls_workbook_path}"
    )

    print(
        f"Elapsed time                    : "
        f"{elapsed_seconds:.2f} s"
    )

    print(
        "=" * 78
    )


# ============================================================
# EXECUTION
# ============================================================

if __name__ == "__main__":

    main()
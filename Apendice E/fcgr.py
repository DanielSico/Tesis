"""
fcgr.py

Codon-FCGR Virus-Host batch analysis.

INPUTS
------

-v / --virus
    Virus FASTA file or directory containing Virus FASTA files.

-h / --host
    Host FASTA file or directory containing Host FASTA files.

-t / --trna
    Plain-text tRNA gene-copy-number table.

-j / --sij
    Plain-text file containing:
        s_GU
        s_UG
        s_IC
        s_IA


OUTPUTS
-------

All outputs are written automatically to:

    results/

Tables:
    unweighted.xlsx
    weighted.xlsx
    SSIM.xlsx
    competition.xlsx

Images:
    results/comparison/
    results/competition/

Only the weighted analysis generates images.


WEIGHTED VIRAL FCGR
-------------------

For sense codon i:

    W_i = sum_j [tGCN_ij * (1 - s_ij)]

    w_i = W_i / max(W)

    inverse weight = 1 / w_i

Weighted viral codon contribution:

    n_i / w_i

The weighted viral FCGR is normalized so that its
64 cells sum to 1.

Host FCGRs remain unweighted relative-frequency matrices.


COMPETITION
-----------

Competition is calculated from the ORIGINAL unweighted
Virus and Host sense-codon relative frequencies:

    competition_i =
        H_i * V_i * (1 / w_i)

The inverse availability term is therefore applied only once.


IMAGES
------

Each Virus-Host comparison generates:

    1 weighted comparison image
    1 competition image

Images are generated at 300 DPI.
"""


# ============================================================
# IMPORTS
# ============================================================

import argparse
import gc
import math
import re
import textwrap
import warnings

from pathlib import Path

import numpy as np
import pandas as pd

from scipy import stats
from statsmodels.stats.multitest import multipletests
from skimage.metrics import structural_similarity as ssim

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from matplotlib.colors import TwoSlopeNorm


# ============================================================
# GENERAL CONSTANTS
# ============================================================

DNA_BASES = "ACGT"

CODON_LENGTH = 3

MATRIX_SIZE = 8

OUTPUT_FOLDER_NAME = "results"

IMAGE_DPI = 300

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
    for codon, amino_acid in GENETIC_CODE.items()
    if amino_acid == "*"
}


SENSE_CODONS = sorted(
    codon
    for codon, amino_acid in GENETIC_CODE.items()
    if amino_acid != "*"
)


SIJ_PARAMETERS = (
    "s_GU",
    "s_UG",
    "s_IC",
    "s_IA",
)


# ============================================================
# FCGR CORNER ASSIGNMENT
# ============================================================

# A = top-left
# C = top-right
# G = bottom-left
# T = bottom-right

CORNER_MAPPING = {
    "A": (0, 0),
    "C": (1, 0),
    "G": (0, 1),
    "T": (1, 1),
}


# ============================================================
# FIXED CODON DISPLAY LAYOUT
# ============================================================

CODON_LAYOUT = (
    "AAA AAC ACA ACC CAA CAC CCA CCC "
    "AAG AAT ACG ACT CAG CAT CCG CCT "
    "AGA AGC ATA ATC CGA CGC CTA CTC "
    "AGG AGT ATG ATT CGG CGT CTG CTT "
    "GAA GAC GCA GCC TAA TAC TCA TCC "
    "GAG GAT GCG GCT TAG TAT TCG TCT "
    "GGA GGC GTA GTC TGA TGC TTA TTC "
    "GGG GGT GTG GTT TGG TGT TTG TTT"
).split()


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
    """
    Return the unique 8x8 FCGR coordinate of one codon.
    """

    x = 0
    y = 0

    for base in codon:

        dx, dy = CORNER_MAPPING[base]

        x = (x << 1) | dx
        y = (y << 1) | dy

    return x, y


CODON_TO_COORDINATE = {
    codon: codon_coordinate(codon)
    for codon in ALL_CODONS
}


COORDINATE_TO_CODON = {
    coordinate: codon
    for codon, coordinate in CODON_TO_COORDINATE.items()
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
    """
    Relative paths are interpreted relative to the directory
    containing this script.
    """

    path = Path(
        path_argument
    )

    if not path.is_absolute():

        path = (
            script_directory()
            / path
        )

    return path.resolve()


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
            f"a file nor a directory:\n"
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
            f"{label} path is not a file:\n"
            f"{path}"
        )


# ============================================================
# FASTA FILE DISCOVERY
# ============================================================

def find_fasta_files(folder):

    fasta_files = sorted(
        (
            path
            for path in folder.iterdir()
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

    current_id = None

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

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if current_id is not None:

                    sequence = (
                        "".join(
                            sequence_parts
                        )
                        .upper()
                    )

                    if not sequence:

                        raise ValueError(
                            f"Empty FASTA sequence for "
                            f"'{current_id}' in:\n"
                            f"{file_path}"
                        )

                    records.append(
                        (
                            current_id,
                            sequence
                        )
                    )

                header = (
                    line[1:]
                    .strip()
                )

                if not header:

                    raise ValueError(
                        f"Empty FASTA header on line "
                        f"{line_number} in:\n"
                        f"{file_path}"
                    )

                current_id = (
                    header
                    .split()[0]
                )

                sequence_parts = []

            else:

                if current_id is None:

                    raise ValueError(
                        f"Sequence data found before "
                        f"the first FASTA header on line "
                        f"{line_number} in:\n"
                        f"{file_path}"
                    )

                sequence_parts.append(
                    line
                )

    if current_id is not None:

        sequence = (
            "".join(
                sequence_parts
            )
            .upper()
        )

        if not sequence:

            raise ValueError(
                f"Empty FASTA sequence for "
                f"'{current_id}' in:\n"
                f"{file_path}"
            )

        records.append(
            (
                current_id,
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

        source_type = "file"

    else:

        fasta_files = (
            find_fasta_files(
                source_path
            )
        )

        source_type = "directory"

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
                sequence_id,
                sequence
            )

        ) in enumerate(
            records,
            start=1
        ):

            if multiple_records:

                sequence_name = (
                    f"{base_name}_"
                    f"{record_number}"
                )

            else:

                sequence_name = (
                    base_name
                )

            invalid_characters = sorted(
                set(sequence)
                - set(DNA_BASES)
            )

            sequences.append(
                {
                    "id":
                        sequence_id,

                    "name":
                        sequence_name,

                    "sequence":
                        sequence,

                    "source_file":
                        file_path,

                    "role":
                        role,

                    "invalid_characters":
                        invalid_characters,
                }
            )

    return {
        "source_path":
            source_path,

        "source_type":
            source_type,

        "fasta_files":
            fasta_files,

        "sequences":
            sequences,
    }


# ============================================================
# CODON EXTRACTION
# ============================================================

def extract_codons(sequence):
    """
    Extract non-overlapping codons beginning at position 0.

    Triplets containing characters other than A, C, G or T
    are skipped.
    """

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
            len(codon) == CODON_LENGTH
            and all(
                base in DNA_BASES
                for base in codon
            )
        ):

            codons.append(
                codon
            )

    return codons


# ============================================================
# tRNA TABLE
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

    text = path.read_text(
        encoding="utf-8-sig"
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
        for codon in ALL_CODONS
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

        if not numeric_count.is_integer():

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
            len(matches),

        "unique_codons":
            len(observed),

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
# sij TABLE
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

    text = path.read_text(
        encoding="utf-8-sig"
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
        for parameter in SIJ_PARAMETERS
        if parameter not in penalties
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
        for base in codon
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
    """
    Calculate:

        W_i =
            sum_j tGCN_ij * (1 - s_ij)

        w_i =
            W_i / max(W_i)

        inverse_i =
            1 / w_i
    """

    edges = (
        build_wobble_edges(
            trna_counts
        )
    )

    raw_weights = {
        codon: 0.0
        for codon in SENSE_CODONS
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
            * (1.0 - sij)
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
# RELATIVE-FREQUENCY FCGR
# ============================================================

def create_relative_frequency_fcgr(
    codons,
    inverse_weights=None
):
    """
    Unweighted:
        contribution = 1

    Weighted Virus:
        contribution = 1 / w_i

    Stop codons:
        contribution = 1

    Matrix is normalized to sum to 1.
    """

    matrix = np.zeros(
        (
            MATRIX_SIZE,
            MATRIX_SIZE
        ),
        dtype=float
    )

    for codon in codons:

        x, y = (
            CODON_TO_COORDINATE[
                codon
            ]
        )

        if inverse_weights is None:

            contribution = 1.0

        elif codon in SENSE_CODONS:

            contribution = float(
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
        np.sum(
            matrix
        )
    )

    if total > 0:

        matrix /= total

    return matrix


# ============================================================
# DIRECT COMPARISON METRICS
# ============================================================

def pairwise_metrics(
    virus_matrix,
    host_matrix
):
    """
    Calculate Virus-Host comparison metrics using standardized
    scientific-library implementations.

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

    spearman_p = float(
        spearman_result.pvalue
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

    pearson_p = float(
        pearson_result.pvalue
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

        "spearman_p":
            spearman_p,

        "pearson_r":
            pearson_r,

        "pearson_p":
            pearson_p,

        "ssim":
            ssim_value,
    }


# ============================================================
# COMPOSITIONAL DISTANCES
# ============================================================

def composition_vector(matrix):

    values = (
        np.asarray(
            matrix,
            dtype=float
        )
        .flatten()
        .copy()
    )

    values = np.nan_to_num(
        values,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    values[
        values < 0
    ] = 0.0

    total = float(
        np.sum(
            values
        )
    )

    if total > 0:

        values /= total

    return values


def distance_metrics(
    virus_matrix,
    host_matrix
):

    p = (
        composition_vector(
            virus_matrix
        )
    )

    q = (
        composition_vector(
            host_matrix
        )
    )

    total_variation = (
        0.5
        * float(
            np.sum(
                np.abs(
                    p - q
                )
            )
        )
    )

    midpoint = (
        0.5
        * (p + q)
    )

    js_divergence = 0.0

    p_mask = (
        p > 0
    )

    if np.any(
        p_mask
    ):

        js_divergence += (
            0.5
            * float(
                np.sum(
                    p[p_mask]
                    * np.log2(
                        p[p_mask]
                        / midpoint[p_mask]
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

        js_divergence += (
            0.5
            * float(
                np.sum(
                    q[q_mask]
                    * np.log2(
                        q[q_mask]
                        / midpoint[q_mask]
                    )
                )
            )
        )

    js_divergence = max(
        js_divergence,
        0.0
    )

    js_distance = float(
        np.sqrt(
            js_divergence
        )
    )

    return {
        "total_variation":
            total_variation,

        "js_distance":
            js_distance,
    }


# ============================================================
# COMPETITION MODEL
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
        np.sum(
            values
        )
    )

    if total > 0:

        values /= total

    return values


def calculate_competition(
    virus_unweighted_matrix,
    host_unweighted_matrix,
    relative_weights,
    weight_matrix,
    sense_mask
):
    """
    Competition uses original unweighted codon frequencies:

        raw_i =
            H_i * V_i / w_i

        directional_i =
            (H_i - V_i) / w_i

    Normalized competition:

        raw_i / (1 / min(w))
    """

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

    theoretical_max = (
        1.0
        / minimum_weight
    )

    scarcity = (
        1.0
        / weight_matrix
    )

    competition_raw = (
        host_frequency
        * virus_frequency
        * scarcity
    )

    competition_normalized = (
        competition_raw
        / theoretical_max
    )

    directional_raw = (
        (
            host_frequency
            - virus_frequency
        )
        * scarcity
    )

    directional_normalized = (
        directional_raw
        / theoretical_max
    )

    directional_normalized = (
        np.clip(
            directional_normalized,
            -1.0,
            1.0
        )
    )

    raw_score = float(
        np.sum(
            competition_raw
        )
    )

    normalized_score = float(
        np.sum(
            competition_normalized
        )
    )

    net_directional_score = float(
        np.sum(
            directional_normalized
        )
    )

    order = (
        np.argsort(
            competition_normalized.ravel()
        )[::-1]
    )

    top_groups = []

    for flat_index in order:

        y, x = (
            np.unravel_index(
                int(
                    flat_index
                ),
                competition_normalized.shape
            )
        )

        value = float(
            competition_normalized[
                y,
                x
            ]
        )

        if value <= 0:
            break

        codon = (
            COORDINATE_TO_CODON[
                (x, y)
            ]
        )

        if codon in STOP_CODONS:
            continue

        top_groups.append(
            {
                "group":
                    codon,

                "competition":
                    value,

                "directional_pressure":
                    float(
                        directional_normalized[
                            y,
                            x
                        ]
                    ),
            }
        )

        if len(top_groups) == 5:
            break

    return {
        "raw_score":
            raw_score,

        "normalized_score":
            normalized_score,

        "net_directional_score":
            net_directional_score,

        "top_groups":
            top_groups,

        "competition_normalized_matrix":
            competition_normalized,

        "directional_normalized_matrix":
            directional_normalized,
    }


# ============================================================
# GC FRACTION
# ============================================================

def calculate_gc_fraction(sequence):

    valid_bases = [
        base
        for base in sequence.upper()
        if base in DNA_BASES
    ]

    if not valid_bases:

        return np.nan

    gc_count = sum(
        1
        for base in valid_bases
        if base in "GC"
    )

    return (
        gc_count
        / len(
            valid_bases
        )
    )


# ============================================================
# SEQUENCE PREPARATION
# ============================================================

def prepare_sequence(
    record,
    inverse_weights=None
):

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

    unweighted_matrix = (
        create_relative_frequency_fcgr(
            codons
        )
    )

    weighted_matrix = None

    if inverse_weights is not None:

        weighted_matrix = (
            create_relative_frequency_fcgr(
                codons,
                inverse_weights=
                    inverse_weights
            )
        )

    return {
        **record,

        "codons":
            codons,

        "unweighted_matrix":
            unweighted_matrix,

        "weighted_matrix":
            weighted_matrix,

        "gc_fraction":
            calculate_gc_fraction(
                record[
                    "sequence"
                ]
            ),
    }


# ============================================================
# RESULT ROW
# ============================================================

def build_result_row(
    virus,
    host,
    virus_matrix,
    host_matrix,
    competition,
    result_label=None,
    result_name=None
):

    direct = (
        pairwise_metrics(
            virus_matrix,
            host_matrix
        )
    )

    distances = (
        distance_metrics(
            virus_matrix,
            host_matrix
        )
    )

    row = {}

    if result_label is not None:

        row[
            "result_label"
        ] = (
            result_label
        )

        row[
            "result_name"
        ] = (
            result_name
        )

    row.update(
        {
            "virus_id":
                virus[
                    "id"
                ],

            "virus_name":
                virus[
                    "name"
                ],

            "host_id":
                host[
                    "id"
                ],

            "host_name":
                host[
                    "name"
                ],

            "vh_spearman_r":
                direct[
                    "spearman_r"
                ],

            "vh_spearman_p":
                direct[
                    "spearman_p"
                ],

            "vh_pearson_r":
                direct[
                    "pearson_r"
                ],

            "vh_pearson_p":
                direct[
                    "pearson_p"
                ],

            "vh_ssim":
                direct[
                    "ssim"
                ],

            "vh_total_variation":
                distances[
                    "total_variation"
                ],

            "vh_js_distance":
                distances[
                    "js_distance"
                ],

            "competition_raw_score":
                competition[
                    "raw_score"
                ],

            "competition_normalized_score":
                competition[
                    "normalized_score"
                ],

            "competition_net_directional_score":
                competition[
                    "net_directional_score"
                ],

            "virus_gc_fraction":
                virus[
                    "gc_fraction"
                ],
        }
    )

    top_groups = (
        competition[
            "top_groups"
        ]
    )

    for rank in range(
        1,
        6
    ):

        if rank <= len(
            top_groups
        ):

            top = (
                top_groups[
                    rank - 1
                ]
            )

            row[
                f"competition_top"
                f"{rank}_group"
            ] = (
                top[
                    "group"
                ]
            )

            row[
                f"competition_top"
                f"{rank}_directional_pressure"
            ] = (
                top[
                    "directional_pressure"
                ]
            )

        else:

            row[
                f"competition_top"
                f"{rank}_group"
            ] = None

            row[
                f"competition_top"
                f"{rank}_directional_pressure"
            ] = np.nan

    return row


# ============================================================
# BATCH ANALYSIS
# ============================================================

def run_batch(
    virus_sequences,
    host_sequences,
    relative_weights,
    inverse_weights
):

    viruses = [
        prepare_sequence(
            record,
            inverse_weights=
                inverse_weights
        )

        for record
        in virus_sequences
    ]

    hosts = [
        prepare_sequence(
            record,
            inverse_weights=None
        )

        for record
        in host_sequences
    ]

    (
        weight_matrix,
        sense_mask

    ) = build_competition_layout(
        relative_weights
    )

    unweighted_rows = []

    weighted_rows = []

    weighted_matrix_records = []

    weighted_image_records = []

    total = (
        len(viruses)
        * len(hosts)
    )

    comparison_number = 0

    for virus in viruses:

        for host in hosts:

            comparison_number += 1

            result_label = (
                f"Comp"
                f"{comparison_number}"
            )

            result_name = (
                f"{virus['name']} "
                f"vs "
                f"{host['name']}"
            )

            print(
                f"[{comparison_number}/"
                f"{total}] "
                f"{result_name}"
            )

            # ------------------------------------------------
            # Competition
            # ------------------------------------------------

            competition = (
                calculate_competition(
                    virus[
                        "unweighted_matrix"
                    ],
                    host[
                        "unweighted_matrix"
                    ],
                    relative_weights,
                    weight_matrix,
                    sense_mask
                )
            )

            # ------------------------------------------------
            # Unweighted result
            # ------------------------------------------------

            unweighted_row = (
                build_result_row(
                    virus=
                        virus,

                    host=
                        host,

                    virus_matrix=
                        virus[
                            "unweighted_matrix"
                        ],

                    host_matrix=
                        host[
                            "unweighted_matrix"
                        ],

                    competition=
                        competition
                )
            )

            # ------------------------------------------------
            # Weighted result
            # ------------------------------------------------

            weighted_row = (
                build_result_row(
                    virus=
                        virus,

                    host=
                        host,

                    virus_matrix=
                        virus[
                            "weighted_matrix"
                        ],

                    host_matrix=
                        host[
                            "unweighted_matrix"
                        ],

                    competition=
                        competition,

                    result_label=
                        result_label,

                    result_name=
                        result_name
                )
            )

            unweighted_rows.append(
                unweighted_row
            )

            weighted_rows.append(
                weighted_row
            )

            # ------------------------------------------------
            # Host - weighted Virus
            # ------------------------------------------------

            comparison_matrix = (
                np.asarray(
                    host[
                        "unweighted_matrix"
                    ],
                    dtype=float
                )
                -
                np.asarray(
                    virus[
                        "weighted_matrix"
                    ],
                    dtype=float
                )
            )

            # ------------------------------------------------
            # Normalized competition matrix
            # ------------------------------------------------

            competition_matrix = (
                np.asarray(
                    competition[
                        "competition_normalized_matrix"
                    ],
                    dtype=float
                )
                .copy()
            )

            weighted_matrix_records.append(
                {
                    "result_label":
                        result_label,

                    "result_name":
                        result_name,

                    "comparison_matrix":
                        comparison_matrix.copy(),

                    "competition_matrix":
                        competition_matrix.copy(),
                }
            )

            weighted_image_records.append(
                {
                    "result_label":
                        result_label,

                    "result_name":
                        result_name,

                    "virus":
                        virus,

                    "host":
                        host,

                    "weighted_virus_matrix":
                        np.asarray(
                            virus[
                                "weighted_matrix"
                            ],
                            dtype=float
                        ).copy(),

                    "host_matrix":
                        np.asarray(
                            host[
                                "unweighted_matrix"
                            ],
                            dtype=float
                        ).copy(),

                    "comparison_matrix":
                        comparison_matrix.copy(),

                    "competition_matrix":
                        competition_matrix.copy(),
                }
            )

    return (
        pd.DataFrame(
            unweighted_rows
        ),

        pd.DataFrame(
            weighted_rows
        ),

        weighted_matrix_records,

        weighted_image_records
    )


# ============================================================
# BENJAMINI-HOCHBERG
# ============================================================

def benjamini_hochberg(values):

    values = np.asarray(
        values,
        dtype=float
    )

    corrected = np.full(
        len(values),
        np.nan,
        dtype=float
    )

    valid = (
        np.isfinite(
            values
        )
        & (values >= 0)
        & (values <= 1)
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


# ============================================================
# FINALIZE RESULT TABLES
# ============================================================

def finalize_results(
    df,
    weighted=False
):

    if df.empty:

        raise ValueError(
            "No Virus-Host comparisons "
            "were generated."
        )

    df = (
        df.copy()
    )

    # --------------------------------------------------------
    # Spearman q
    # --------------------------------------------------------

    df[
        "vh_spearman_q"
    ] = benjamini_hochberg(
        pd.to_numeric(
            df[
                "vh_spearman_p"
            ],
            errors="coerce"
        )
        .to_numpy(
            dtype=float
        )
    )

    # --------------------------------------------------------
    # Pearson q
    # --------------------------------------------------------

    df[
        "vh_pearson_q"
    ] = benjamini_hochberg(
        pd.to_numeric(
            df[
                "vh_pearson_p"
            ],
            errors="coerce"
        )
        .to_numpy(
            dtype=float
        )
    )

    # --------------------------------------------------------
    # Competition rank
    # --------------------------------------------------------

    scores = (
        pd.to_numeric(
            df[
                "competition_normalized_score"
            ],
            errors="coerce"
        )
    )

    ranks = (
        scores.rank(
            method="min",
            ascending=False,
            na_option="bottom"
        )
    )

    df[
        "competition_rank"
    ] = (
        ranks.astype(
            "Int64"
        )
    )

    core_columns = [
        "virus_id",
        "virus_name",
        "host_id",
        "host_name",

        "vh_spearman_r",
        "vh_spearman_p",
        "vh_spearman_q",

        "vh_pearson_r",
        "vh_pearson_p",
        "vh_pearson_q",

        "vh_ssim",
        "vh_total_variation",
        "vh_js_distance",

        "competition_raw_score",
        "competition_normalized_score",
        "competition_rank",
        "competition_net_directional_score",

        "competition_top1_group",
        "competition_top1_directional_pressure",

        "competition_top2_group",
        "competition_top2_directional_pressure",

        "competition_top3_group",
        "competition_top3_directional_pressure",

        "competition_top4_group",
        "competition_top4_directional_pressure",

        "competition_top5_group",
        "competition_top5_directional_pressure",

        "virus_gc_fraction",
    ]

    if weighted:

        columns = [
            "result_label",
            "result_name",
            *core_columns,
        ]

    else:

        columns = (
            core_columns
        )

    return (
        df[
            columns
        ]
        .copy()
    )


# ============================================================
# PAIRWISE MATRIX CORRELATIONS
# ============================================================

def matrix_pearson_r(
    matrix_a,
    matrix_b
):

    a = np.asarray(
        matrix_a,
        dtype=float
    ).flatten()

    b = np.asarray(
        matrix_b,
        dtype=float
    ).flatten()

    with warnings.catch_warnings():

        warnings.simplefilter(
            "ignore"
        )

        result = (
            stats.pearsonr(
                a,
                b
            )
        )

    return float(
        result.statistic
    )


def build_pairwise_matrix(
    matrix_records,
    result_df,
    matrix_key,
    order_column
):

    if (
        len(matrix_records)
        != len(result_df)
    ):

        raise ValueError(
            "Matrix-record count does not match "
            "result-table row count."
        )

    order_values = (
        pd.to_numeric(
            result_df[
                order_column
            ],
            errors="coerce"
        )
        .to_numpy(
            dtype=float
        )
    )

    def sort_key(index):

        value = (
            order_values[
                index
            ]
        )

        if np.isfinite(
            value
        ):

            return (
                0,
                -float(
                    value
                ),
                index
            )

        return (
            1,
            0.0,
            index
        )

    order = sorted(
        range(
            len(
                matrix_records
            )
        ),
        key=sort_key
    )

    ordered = [
        matrix_records[
            index
        ]
        for index
        in order
    ]

    labels = [
        (
            f"{record['result_label']}: "
            f"{record['result_name']}"
        )
        for record
        in ordered
    ]

    number_of_matrices = (
        len(
            ordered
        )
    )

    values = np.full(
        (
            number_of_matrices,
            number_of_matrices
        ),
        np.nan,
        dtype=float
    )

    for row_index in range(
        number_of_matrices
    ):

        values[
            row_index,
            row_index
        ] = 1.0

        for col_index in range(
            row_index + 1,
            number_of_matrices
        ):

            value = (
                matrix_pearson_r(
                    ordered[
                        row_index
                    ][
                        matrix_key
                    ],
                    ordered[
                        col_index
                    ][
                        matrix_key
                    ]
                )
            )

            values[
                row_index,
                col_index
            ] = value

            values[
                col_index,
                row_index
            ] = value

    return pd.DataFrame(
        values,
        index=labels,
        columns=labels
    )


# ============================================================
# IMAGE HELPERS
# ============================================================

def safe_filename(text):

    text = str(
        text
    )

    text = re.sub(
        r'[<>:"/\\|?*\x00-\x1F]+',
        "_",
        text
    )

    text = re.sub(
        r"\s+",
        " ",
        text
    ).strip()

    return text[:220]


def codon_layout_matrix_values(matrix):
    """
    Reorder matrix values according to the explicit
    display layout.
    """

    values = np.array(
        [
            float(
                matrix[
                    CODON_TO_COORDINATE[
                        codon
                    ][1],
                    CODON_TO_COORDINATE[
                        codon
                    ][0]
                ]
            )

            for codon
            in CODON_LAYOUT
        ],
        dtype=float
    ).reshape(
        MATRIX_SIZE,
        MATRIX_SIZE
    )

    return values


# ============================================================
# IMAGE RENDERING
# ============================================================

def render_weighted_fcgr_figure(
    virus,
    host,
    weighted_virus_matrix,
    host_matrix,
    comparison_matrix,
    competition_matrix,
    output_path,
    mode
):
    """
    Generate one three-panel FCGR PNG.

    Memory management:
        - figure explicitly cleared
        - figure explicitly closed
        - all Matplotlib figures closed
        - temporary arrays deleted
        - gc.collect() forced
    """

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    virus_values = (
        codon_layout_matrix_values(
            weighted_virus_matrix
        )
    )

    host_values = (
        codon_layout_matrix_values(
            host_matrix
        )
    )

    if mode == "comparison":

        third_values = (
            codon_layout_matrix_values(
                comparison_matrix
            )
        )

        third_title = (
            "difference (Host - Virus)"
        )

    elif mode == "competition":

        third_values = (
            codon_layout_matrix_values(
                competition_matrix
            )
        )

        third_title = (
            "competition"
        )

    else:

        raise ValueError(
            f"Unknown image mode: "
            f"{mode}"
        )

    titles = [
        virus[
            "name"
        ],

        host[
            "name"
        ],

        third_title,
    ]

    fig = None

    try:

        fig, axes = plt.subplots(
            1,
            3,
            figsize=(
                20,
                7.2
            ),
            layout="constrained"
        )

        max_frequency = max(
            float(
                np.max(
                    virus_values
                )
            ),
            float(
                np.max(
                    host_values
                )
            ),
            1e-12
        )

        matrices = [
            virus_values,
            host_values,
            third_values,
        ]

        for (
            panel_index,
            (
                ax,
                values,
                title
            )

        ) in enumerate(
            zip(
                axes,
                matrices,
                titles
            )
        ):

            # ------------------------------------------------
            # Difference panel
            # ------------------------------------------------

            if (
                panel_index == 2
                and mode == "comparison"
            ):

                limit = max(
                    0.08,
                    float(
                        np.abs(
                            values
                        ).max()
                    )
                )

                image = ax.imshow(
                    values,
                    cmap="RdBu",
                    norm=TwoSlopeNorm(
                        vmin=-limit,
                        vcenter=0.0,
                        vmax=limit
                    ),
                    interpolation="nearest"
                )

            # ------------------------------------------------
            # Virus / Host / Competition
            # ------------------------------------------------

            else:

                if panel_index < 2:

                    vmax = (
                        max_frequency
                    )

                else:

                    vmax = max(
                        float(
                            np.max(
                                values
                            )
                        ),
                        1e-12
                    )

                image = ax.imshow(
                    values,
                    cmap="Blues",
                    vmin=0.0,
                    vmax=vmax,
                    interpolation="nearest"
                )

            # ------------------------------------------------
            # Title
            # ------------------------------------------------

            ax.set_title(
                textwrap.fill(
                    title,
                    width=48,
                    break_long_words=False,
                    break_on_hyphens=False
                ),
                fontsize=13,
                pad=14
            )

            ax.set_xticks(
                []
            )

            ax.set_yticks(
                []
            )

            # ------------------------------------------------
            # Cell labels
            # ------------------------------------------------

            for (
                index,
                codon

            ) in enumerate(
                CODON_LAYOUT
            ):

                row = (
                    index // 8
                )

                column = (
                    index % 8
                )

                value = float(
                    values[
                        row,
                        column
                    ]
                )

                if (
                    panel_index == 2
                    and mode == "competition"
                ):

                    shown = (
                        f"{value * 1e4:.2f}"
                    )

                else:

                    shown = (
                        f"{value:.3f}"
                    )

                if shown == "-0.000":

                    shown = (
                        "0.000"
                    )

                rgba = (
                    image.cmap(
                        image.norm(
                            value
                        )
                    )
                )

                luminance = (
                    0.2126 * rgba[0]
                    + 0.7152 * rgba[1]
                    + 0.0722 * rgba[2]
                )

                if luminance < 0.48:

                    text_color = (
                        "white"
                    )

                else:

                    text_color = (
                        "black"
                    )

                ax.text(
                    column,
                    row,
                    f"{codon}\n"
                    f"{shown}",
                    ha="center",
                    va="center",
                    fontsize=9,
                    color=text_color
                )

            # ------------------------------------------------
            # Colorbar
            # ------------------------------------------------

            colorbar = (
                fig.colorbar(
                    image,
                    ax=ax,
                    shrink=0.78,
                    pad=0.025
                )
            )

            colorbar.ax.tick_params(
                labelsize=9
            )

            if panel_index < 2:

                colorbar.set_label(
                    "Relative frequency",
                    fontsize=10
                )

            elif mode == "comparison":

                colorbar.set_label(
                    "Host - Virus",
                    fontsize=10
                )

            else:

                colorbar.set_label(
                    "Normalized contribution",
                    fontsize=10
                )

        # ----------------------------------------------------
        # Figure note
        # ----------------------------------------------------

        if mode == "comparison":

            note = (
                "Virus: inverse tRNA-availability weighted "
                "relative frequency. "
                "Host: unweighted relative frequency. "
                "Difference: Host - weighted Virus."
            )

        else:

            note = (
                "Virus: inverse tRNA-availability weighted "
                "relative frequency. "
                "Host: unweighted relative frequency. "
                "Competition cell labels × 10,000; "
                "colorbar shown in original normalized units."
            )

        fig.supxlabel(
            note,
            fontsize=10
        )

        # ----------------------------------------------------
        # Save
        # ----------------------------------------------------

        fig.savefig(
            output_path,
            dpi=IMAGE_DPI
        )

    finally:

        # ====================================================
        # EXPLICIT MEMORY CLEANUP
        # ====================================================

        if fig is not None:

            fig.clear()

            plt.close(
                fig
            )

        plt.close(
            "all"
        )

        del virus_values
        del host_values
        del third_values

        gc.collect()


# ============================================================
# EXPORT WEIGHTED IMAGES
# ============================================================

def export_weighted_fcgr_images(
    weighted_image_records,
    output_directory
):
    """
    Every execution regenerates all images.

    Existing files are overwritten.

    No resume behavior is used.
    """

    comparison_directory = (
        output_directory
        / "comparison"
    )

    competition_directory = (
        output_directory
        / "competition"
    )

    comparison_directory.mkdir(
        parents=True,
        exist_ok=True
    )

    competition_directory.mkdir(
        parents=True,
        exist_ok=True
    )

    total = len(
        weighted_image_records
    )

    print()

    print(
        "=" * 78
    )

    print(
        "EXPORTING WEIGHTED FCGR IMAGES"
    )

    print(
        "=" * 78
    )

    print(
        f"Comparison figures              : "
        f"{total}"
    )

    print(
        f"Competition figures             : "
        f"{total}"
    )

    print(
        f"Total PNG images                : "
        f"{total * 2}"
    )

    print(
        f"Image resolution                : "
        f"{IMAGE_DPI} DPI"
    )

    print()

    for (
        index,
        record

    ) in enumerate(
        weighted_image_records,
        start=1
    ):

        result_label = (
            record[
                "result_label"
            ]
        )

        result_name = (
            record[
                "result_name"
            ]
        )

        filename_stem = (
            safe_filename(
                f"{result_label} -- "
                f"{result_name}"
            )
        )

        comparison_path = (
            comparison_directory
            / (
                f"{filename_stem}"
                f".png"
            )
        )

        competition_path = (
            competition_directory
            / (
                f"{filename_stem}"
                f".png"
            )
        )

        # ----------------------------------------------------
        # Comparison image
        # ----------------------------------------------------

        render_weighted_fcgr_figure(
            virus=
                record[
                    "virus"
                ],

            host=
                record[
                    "host"
                ],

            weighted_virus_matrix=
                record[
                    "weighted_virus_matrix"
                ],

            host_matrix=
                record[
                    "host_matrix"
                ],

            comparison_matrix=
                record[
                    "comparison_matrix"
                ],

            competition_matrix=
                record[
                    "competition_matrix"
                ],

            output_path=
                comparison_path,

            mode=
                "comparison"
        )

        # ----------------------------------------------------
        # Explicit cleanup between the two images
        # ----------------------------------------------------

        plt.close(
            "all"
        )

        gc.collect()

        # ----------------------------------------------------
        # Competition image
        # ----------------------------------------------------

        render_weighted_fcgr_figure(
            virus=
                record[
                    "virus"
                ],

            host=
                record[
                    "host"
                ],

            weighted_virus_matrix=
                record[
                    "weighted_virus_matrix"
                ],

            host_matrix=
                record[
                    "host_matrix"
                ],

            comparison_matrix=
                record[
                    "comparison_matrix"
                ],

            competition_matrix=
                record[
                    "competition_matrix"
                ],

            output_path=
                competition_path,

            mode=
                "competition"
        )

        # ----------------------------------------------------
        # Explicit cleanup before next pair
        # ----------------------------------------------------

        plt.close(
            "all"
        )

        gc.collect()

        print(
            f"[{index}/{total}] "
            f"{result_label}: "
            f"{result_name}"
        )


# ============================================================
# CONSOLE SEQUENCE SUMMARY
# ============================================================

def print_sequence_summary(
    title,
    result
):

    print()

    print(
        "=" * 78
    )

    print(
        title
    )

    print(
        "=" * 78
    )

    print(
        f"Input                           : "
        f"{result['source_path']}"
    )

    print(
        f"Input type                      : "
        f"{result['source_type']}"
    )

    print(
        f"FASTA files                     : "
        f"{len(result['fasta_files']):,}"
    )

    print(
        f"Sequences loaded                : "
        f"{len(result['sequences']):,}"
    )

    print(
        "-" * 78
    )

    for (
        index,
        record

    ) in enumerate(
        result[
            "sequences"
        ],
        start=1
    ):

        print(
            f"{index:>3}. "
            f"{record['name']}"
        )

        print(
            f"     ID: "
            f"{record['id']}"
        )

        print(
            f"     Length: "
            f"{len(record['sequence']):,} nt"
        )

        if record[
            "invalid_characters"
        ]:

            print(
                "     WARNING: non-ACGT characters: "
                + ", ".join(
                    record[
                        "invalid_characters"
                    ]
                )
            )


# ============================================================
# MAIN
# ============================================================

def main():

    # --------------------------------------------------------
    # -h belongs to Host.
    # Help is therefore --help.
    # --------------------------------------------------------

    parser = argparse.ArgumentParser(
        description=(
            "Unweighted and inverse-tRNA-availability "
            "weighted codon-FCGR Virus-Host analysis."
        ),
        add_help=False
    )

    # --------------------------------------------------------
    # Help
    # --------------------------------------------------------

    parser.add_argument(
        "--help",
        action="help",
        help=(
            "Show this help message and exit."
        )
    )

    # --------------------------------------------------------
    # Virus
    # --------------------------------------------------------

    parser.add_argument(
        "-v",
        "--virus",
        required=True,
        help=(
            "Virus FASTA file or directory containing "
            "Virus FASTA files."
        )
    )

    # --------------------------------------------------------
    # Host
    # --------------------------------------------------------

    parser.add_argument(
        "-h",
        "--host",
        required=True,
        help=(
            "Host FASTA file or directory containing "
            "Host FASTA files."
        )
    )

    # --------------------------------------------------------
    # tRNA
    # --------------------------------------------------------

    parser.add_argument(
        "-t",
        "--trna",
        required=True,
        help=(
            "Plain-text tRNA codon/copy-number file."
        )
    )

    # --------------------------------------------------------
    # sij
    # --------------------------------------------------------

    parser.add_argument(
        "-j",
        "--sij",
        required=True,
        help=(
            "Plain-text file containing "
            "s_GU, s_UG, s_IC and s_IA."
        )
    )

    args = (
        parser.parse_args()
    )

    # ========================================================
    # INPUT PATHS
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

    # ========================================================
    # OUTPUT DIRECTORY
    # ========================================================

    output_directory = (
        script_directory()
        / OUTPUT_FOLDER_NAME
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True
    )

    # ========================================================
    # VALIDATE INPUTS
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

    # ========================================================
    # LOAD INPUTS
    # ========================================================

    virus_result = (
        load_sequence_source(
            virus_source,
            role="virus"
        )
    )

    host_result = (
        load_sequence_source(
            host_source,
            role="host"
        )
    )

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

    # ========================================================
    # INPUT SUMMARY
    # ========================================================

    print_sequence_summary(
        "VIRUS SEQUENCES",
        virus_result
    )

    print_sequence_summary(
        "HOST SEQUENCES",
        host_result
    )

    print()

    print(
        "=" * 78
    )

    print(
        "tRNA AND sij MODEL"
    )

    print(
        "=" * 78
    )

    print(
        f"tRNA file                       : "
        f"{trna_path}"
    )

    print(
        f"Codon/count pairs read          : "
        f"{trna_result['pairs_read']}"
    )

    print(
        f"Unique codons recognized        : "
        f"{trna_result['unique_codons']}/64"
    )

    print(
        f"Total tRNA copies               : "
        f"{trna_result['total_copies']}"
    )

    if trna_result[
        "duplicates"
    ]:

        print(
            "Duplicate codon entries          : "
            + ", ".join(
                trna_result[
                    "duplicates"
                ]
            )
            + " (summed)"
        )

    print(
        f"sij file                        : "
        f"{sij_path}"
    )

    for parameter in (
        SIJ_PARAMETERS
    ):

        print(
            f"{parameter:<32}: "
            f"{sij_values[parameter]:.6f}"
        )

    # ========================================================
    # tRNA AVAILABILITY
    # ========================================================

    (
        raw_weights,
        relative_weights,
        inverse_weights,
        edges

    ) = calculate_relative_codon_weights(
        trna_result[
            "counts"
        ],
        sij_values
    )

    print()

    print(
        "=" * 78
    )

    print(
        "INVERSE tRNA-AVAILABILITY WEIGHTING"
    )

    print(
        "=" * 78
    )

    print(
        f"Positive-copy decoding edges    : "
        f"{len(edges)}"
    )

    print(
        f"Maximum raw availability        : "
        f"{max(raw_weights.values()):.6f}"
    )

    print(
        f"Minimum relative availability   : "
        f"{min(relative_weights.values()):.6f}"
    )

    print(
        f"Maximum inverse penalty         : "
        f"{max(inverse_weights.values()):.6f}"
    )

    print(
        "Weighted Virus contribution     : "
        "1 / w_i"
    )

    print(
        "Stop-codon contribution         : "
        "1.0"
    )

    # ========================================================
    # BATCH SUMMARY
    # ========================================================

    number_of_viruses = len(
        virus_result[
            "sequences"
        ]
    )

    number_of_hosts = len(
        host_result[
            "sequences"
        ]
    )

    comparison_count = (
        number_of_viruses
        * number_of_hosts
    )

    print()

    print(
        "=" * 78
    )

    print(
        "CODON-FCGR BATCH"
    )

    print(
        "=" * 78
    )

    print(
        f"Codon length                    : "
        f"{CODON_LENGTH} nucleotides"
    )

    print(
        f"FCGR matrix                     : "
        f"{MATRIX_SIZE} x {MATRIX_SIZE}"
    )

    print(
        "Frequency representation        : "
        "relative frequency"
    )

    print(
        f"Virus sequences                 : "
        f"{number_of_viruses}"
    )

    print(
        f"Host sequences                  : "
        f"{number_of_hosts}"
    )

    print(
        f"Virus x Host comparisons        : "
        f"{comparison_count}"
    )

    print(
        f"Output directory                : "
        f"{output_directory}"
    )

    print()

    # ========================================================
    # RUN ANALYSIS
    # ========================================================

    (
        unweighted_results,
        weighted_results,
        weighted_matrix_records,
        weighted_image_records

    ) = run_batch(
        virus_result[
            "sequences"
        ],
        host_result[
            "sequences"
        ],
        relative_weights,
        inverse_weights
    )

    # ========================================================
    # FINALIZE TABLES
    # ========================================================

    unweighted_results = (
        finalize_results(
            unweighted_results,
            weighted=False
        )
    )

    weighted_results = (
        finalize_results(
            weighted_results,
            weighted=True
        )
    )

    # ========================================================
    # SSIM-ORDERED MATRIX
    # ========================================================

    ssim_matrix = (
        build_pairwise_matrix(
            weighted_matrix_records,
            weighted_results,
            matrix_key=
                "comparison_matrix",
            order_column=
                "vh_ssim"
        )
    )

    # ========================================================
    # COMPETITION-ORDERED MATRIX
    # ========================================================

    competition_matrix = (
        build_pairwise_matrix(
            weighted_matrix_records,
            weighted_results,
            matrix_key=
                "competition_matrix",
            order_column=
                "competition_normalized_score"
        )
    )

    # ========================================================
    # OUTPUT PATHS
    # ========================================================

    unweighted_path = (
        output_directory
        / "unweighted.xlsx"
    )

    weighted_path = (
        output_directory
        / "weighted.xlsx"
    )

    ssim_path = (
        output_directory
        / "SSIM.xlsx"
    )

    competition_path = (
        output_directory
        / "competition.xlsx"
    )

    # ========================================================
    # WRITE EXCEL FILES
    # ========================================================

    unweighted_results.to_excel(
        unweighted_path,
        sheet_name="Unweighted",
        index=False,
        freeze_panes=(
            1,
            0
        )
    )

    weighted_results.to_excel(
        weighted_path,
        sheet_name="Weighted",
        index=False,
        freeze_panes=(
            1,
            0
        )
    )

    ssim_matrix.to_excel(
        ssim_path,
        sheet_name="SSIM Pearson r",
        index=True,
        index_label="Result",
        freeze_panes=(
            1,
            1
        )
    )

    competition_matrix.to_excel(
        competition_path,
        sheet_name="Competition Pearson r",
        index=True,
        index_label="Result",
        freeze_panes=(
            1,
            1
        )
    )

    # ========================================================
    # CLEAN BEFORE IMAGE GENERATION
    # ========================================================

    plt.close(
        "all"
    )

    gc.collect()

    # ========================================================
    # EXPORT WEIGHTED IMAGES ONLY
    # ========================================================

    export_weighted_fcgr_images(
        weighted_image_records=
            weighted_image_records,

        output_directory=
            output_directory
    )

    # ========================================================
    # FINAL MEMORY CLEANUP
    # ========================================================

    plt.close(
        "all"
    )

    gc.collect()

    # ========================================================
    # FINAL REPORT
    # ========================================================

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
        f"Comparisons                     : "
        f"{len(weighted_results)}"
    )

    print(
        f"unweighted.xlsx                 : "
        f"{unweighted_path}"
    )

    print(
        f"weighted.xlsx                   : "
        f"{weighted_path}"
    )

    print(
        f"SSIM.xlsx                       : "
        f"{ssim_path}"
    )

    print(
        f"competition.xlsx                : "
        f"{competition_path}"
    )

    print(
        f"Comparison PNGs                 : "
        f"{output_directory / 'comparison'}"
    )

    print(
        f"Competition PNGs                : "
        f"{output_directory / 'competition'}"
    )

    print(
        f"Image resolution                : "
        f"{IMAGE_DPI} DPI"
    )

    print(
        "=" * 78
    )


# ============================================================
# EXECUTION
# ============================================================

if __name__ == "__main__":

    main()
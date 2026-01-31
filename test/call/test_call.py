import os
import pytest

from ..helpers.integration import run_integration_test


COMMAND = os.path.basename(os.path.dirname(__file__))
GZIP = ".gz"
OVERWRITE_OUTPUT = False
INFILE_BAM = f"test/{COMMAND}/input/HG002_chr1_MATERNAL_1-500000.bam"
INFILE_FA = f"test/{COMMAND}/input/HG002_chr1_MATERNAL_1-500000.fa.gz"
INFILE_BED = f"test/{COMMAND}/input/HG002_chr1_MATERNAL.bed"
INFILE_CFG = f"test/{COMMAND}/input/config.toml"
PILEUP_TYPES = (
    "cov",
    "mismatch",
    "mapq",
    "insertion",
    "deletion",
    "softclip",
    "ident",
)


@pytest.mark.parametrize(
    ["expected", "added_flags"],
    [
        # Default w/config
        (
            f"test/{COMMAND}/expected/calls_default.bed{GZIP}",
            [],
        ),
        # Fasta
        (
            f"test/{COMMAND}/expected/calls_with_fasta.bed{GZIP}",
            ["-f", INFILE_FA],
        ),
        # Plot
        (
            # Doesn't check image.
            f"test/{COMMAND}/expected/calls_plot.bed{GZIP}",
            ["-d", f"test/{COMMAND}/expected/plot"],
        ),
        # Plot + overlap_calls
        (
            f"test/{COMMAND}/expected/calls_plot_overlap_calls.bed{GZIP}",
            ["-d", f"test/{COMMAND}/expected/plot_overlap_calls", "--overlap_calls"],
        ),
        # Plot + ignore_mtypes
        (
            f"test/{COMMAND}/expected/calls_ignore_mtypes.bed{GZIP}",
            [
                "--ignore_mtypes",
                "dinucleotide",
                "-f",
                INFILE_FA,
                "-d",
                f"test/{COMMAND}/expected/plot_ignore_mtypes",
                "--overlap_calls",
            ],
        ),
        # Plot + builtin_tracks
        (
            f"test/{COMMAND}/expected/calls_plot_builtin_tracks.bed{GZIP}",
            [
                "-d",
                f"test/{COMMAND}/expected/plot_builtin_tracks",
                "-f",
                INFILE_FA,
                "--add_builtin_tracks",
                "ident",
                "mapq",
            ],
        ),
        # Plot + builtin_tracks + ident_breakpoints
        (
            f"test/{COMMAND}/expected/calls_plot_builtin_tracks_ident_breakpoints.bed{GZIP}",
            [
                "-d",
                f"test/{COMMAND}/expected/plot_builtin_tracks_ident_breakpoints",
                "-f",
                INFILE_FA,
                "--add_builtin_tracks",
                "ident",
                "--ident_breakpoints",
                f"test/{COMMAND}/input/breakpoints.tsv",
            ],
        ),
        # Plot + tracks
        (
            f"test/{COMMAND}/expected/calls_plot_tracks.bed{GZIP}",
            [
                "-d",
                f"test/{COMMAND}/expected/plot_tracks",
                "--tracks",
                f"test/{COMMAND}/input/tracks.bed",
            ],
        ),
        # Pileup wigs
        (
            f"test/{COMMAND}/expected/calls_pileup.bed{GZIP}",
            [
                "--output_pileup_dir",
                f"test/{COMMAND}/expected/pileup",
                "--add_pileup_data",
                *PILEUP_TYPES,
            ],
        ),
        # Pileup bigWigs
        (
            f"test/{COMMAND}/expected/calls_pileup_bigwig.bed{GZIP}",
            [
                "--output_pileup_dir",
                f"test/{COMMAND}/expected/pileup",
                "-f",
                INFILE_FA,
                "--add_pileup_data",
                *PILEUP_TYPES,
            ],
        ),
    ],
)
def test_call(expected, added_flags):
    cmd = [
        "python",
        "-m",
        "nucflag.main",
        COMMAND,
        "-i",
        INFILE_BAM,
        "-b",
        INFILE_BED,
        "-c",
        INFILE_CFG,
        *added_flags,
    ]
    run_integration_test(
        *cmd, expected_output=expected, overwrite_output=OVERWRITE_OUTPUT
    )

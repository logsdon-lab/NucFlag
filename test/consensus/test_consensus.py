import os
import pytest

from ..helpers.integration import run_integration_test


COMMAND = os.path.basename(os.path.dirname(__file__))
GZIP = ".gz"
OVERWRITE_OUTPUT = False
INFILES = (
    f"test/{COMMAND}/input/HG002_HMM-Flagger_v1.1.0_HiFi.bed.gz",
    f"test/{COMMAND}/input/HG002_NucFlag_v0.3.7_HiFi.bed.gz",
    f"test/{COMMAND}/input/HG002_NucFlag_v1.0_ONT.bed.gz",
)


@pytest.mark.parametrize(
    ["infiles", "expected", "added_flags"],
    [
        # Test any overlap
        (INFILES, f"test/{COMMAND}/expected/any_overlap.bed{GZIP}", []),
        # Number of overlap
        # Percent overlap
        (INFILES, f"test/{COMMAND}/expected/perc_overlap.bed{GZIP}", ["-p", "0.5"]),
        # Test overlap A, B, A_B
        *[
            (
                INFILES,
                f"test/{COMMAND}/expected/type_{typ}_overlap.bed{GZIP}",
                ["-t", typ, "-p", "0.5"],
            )
            for typ in ("A", "B", "A_B")
        ],
        # Filter calls
        (INFILES, f"test/{COMMAND}/expected/filter_calls.bed{GZIP}", ["-f", "Dup"]),
        # Relative to call
        *[
            (
                INFILES,
                f"test/{COMMAND}/expected/relative_to_{i}.bed{GZIP}",
                ["-r", str(i)],
            )
            for i in range(len(INFILES))
        ],
    ],
)
def test_consensus(infiles, expected, added_flags):
    cmd = ["python", "-m", "nucflag.main", COMMAND, "-i", *infiles, *added_flags]
    run_integration_test(
        *cmd, expected_output=expected, overwrite_output=OVERWRITE_OUTPUT
    )

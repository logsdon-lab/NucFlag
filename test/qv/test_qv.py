import os
import pytest

from ..helpers.integration import run_integration_test


COMMAND = os.path.basename(os.path.dirname(__file__))
GZIP = ".gz"
OVERWRITE_OUTPUT = False


@pytest.mark.parametrize(
    ["infile", "expected", "added_flags"],
    [
        (
            "test/ideogram/input/HG002.bed.gz",
            f"test/{COMMAND}/expected/default.bed{GZIP}",
            [],
        ),
        (
            "test/ideogram/input/HG002.bed.gz",
            f"test/{COMMAND}/expected/exclude_homopolymers.bed{GZIP}",
            ["-c", "homopolymer"],
        ),
    ],
)
def test_qv(infile, expected, added_flags):
    cmd = ["python", "-m", "nucflag.main", COMMAND, "-i", infile, *added_flags]
    run_integration_test(
        *cmd, expected_output=expected, overwrite_output=OVERWRITE_OUTPUT
    )

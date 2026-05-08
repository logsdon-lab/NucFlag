import os
import pytest

from ..helpers.integration import run_integration_test


COMMAND = os.path.basename(os.path.dirname(__file__))
GZIP = ""
OVERWRITE_OUTPUT = False


@pytest.mark.parametrize(
    ["infile", "bed", "expected", "added_flags"],
    [
        (
            "test/consensus/input/HG002_NucFlag_v1.0_ONT.bed.gz",
            f"test/{COMMAND}/input/HG002_cen.bed",
            f"test/{COMMAND}/expected/status_name.tsv{GZIP}",
            ["-g", "name"],
        ),
        (
            "test/consensus/input/HG002_NucFlag_v1.0_ONT.bed.gz",
            f"test/{COMMAND}/input/HG002_cen_region.bed",
            f"test/{COMMAND}/expected/status_region.bed{GZIP}",
            [],
        ),
        (
            "test/consensus/input/HG002_NucFlag_v1.0_ONT.bed.gz",
            f"test/{COMMAND}/input/HG002_cen_region.bed",
            f"test/{COMMAND}/expected/status_region_count.bed{GZIP}",
            ["-m", "count"],
        ),
    ],
)
def test_status(infile, bed, expected, added_flags):
    cmd = [
        "python",
        "-m",
        "nucflag.main",
        COMMAND,
        "-i",
        infile,
        "-b",
        bed,
        *added_flags,
    ]
    run_integration_test(
        *cmd, expected_output=expected, overwrite_output=OVERWRITE_OUTPUT
    )

import os
import pytest

from ..helpers.integration import run_integration_test


COMMAND = os.path.basename(os.path.dirname(__file__))
GZIP = ""
OVERWRITE_OUTPUT = False


@pytest.mark.parametrize(
    ["infile", "bed", "expected", "added_flags"],
    [
        # Status by group name in bed
        (
            "test/consensus/input/HG002_NucFlag_v1.0_ONT.bed.gz",
            f"test/{COMMAND}/input/HG002_cen.bed",
            f"test/{COMMAND}/expected/status_name.tsv{GZIP}",
            ["-g", "name"],
        ),
        # Status by region in bed
        (
            "test/consensus/input/HG002_NucFlag_v1.0_ONT.bed.gz",
            f"test/{COMMAND}/input/HG002_cen_region.bed",
            f"test/{COMMAND}/expected/status_region.bed{GZIP}",
            [],
        ),
        # Status by region with count in bed
        (
            "test/consensus/input/HG002_NucFlag_v1.0_ONT.bed.gz",
            f"test/{COMMAND}/input/HG002_cen_region.bed",
            f"test/{COMMAND}/expected/status_region_count.bed{GZIP}",
            ["-m", "count"],
        ),
        # Status by region with overlapping regions in bed
        (
            f"test/{COMMAND}/input/CHM13_chr4_calls.bed.gz",
            f"test/{COMMAND}/input/CHM13_chr4_alr_ovl.bed",
            f"test/{COMMAND}/expected/status_ovl_regions.bed{GZIP}",
            [],
        ),
        # Genome-wide status with no bed
        (
            f"test/{COMMAND}/input/HG00097_hap1_NucFlag_v1.0.0-a2_hifi.bed.gz",
            None,
            f"test/{COMMAND}/expected/status_all.tsv{GZIP}",
            ["-g", "all"],
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
        *added_flags,
    ]
    if bed:
        cmd.extend(["-b", bed])

    run_integration_test(
        *cmd, expected_output=expected, overwrite_output=OVERWRITE_OUTPUT
    )

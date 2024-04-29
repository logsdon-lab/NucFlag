import pytest
import subprocess


@pytest.mark.parametrize(
    ["bam", "bed", "expected", "config"],
    [
        (
            "test/HG00096_hifi_test.bam",
            "test/test.bed",
            "test/expected.bed",
            "test/config.toml",
        )
    ],
)
def test_identify_misassemblies(bam: str, bed: str, expected: str, config: str):
    process = subprocess.run(
        ["nucflag", "-i", bam, "-b", bed, "-c", config],
        capture_output=True,
        check=True,
    )
    res = [line.split("\t") for line in process.stdout.decode().split("\n") if line]
    with open(expected, "rt") as exp_res_fh:
        exp_res = [line.strip().split("\t") for line in exp_res_fh.readlines() if line]
        assert res == exp_res

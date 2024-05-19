import pytest
import subprocess


@pytest.mark.parametrize(
    ["bam", "bed", "expected", "config"],
    [
        (
            "test/standard/HG00096_hifi.bam",
            "test/standard/region.bed",
            "test/standard/expected.bed",
            tuple(["-c", "test/config.toml"]),
        ),
        (
            "test/ignored/HG00731_hifi.bam",
            "test/ignored/region.bed",
            "test/ignored/expected.bed",
            tuple(["-c", "test/config.toml", "--ignore_regions", "test/ignored/ignore.bed"]),
        )
    ],
)
def test_identify_misassemblies(bam: str, bed: str, expected: str, config: tuple[str]):
    process = subprocess.run(
        ["python", "-m", "nucflag.main", "-i", bam, "-b", bed, *config],
        capture_output=True,
        check=True,
    )
    res = [line.split("\t") for line in process.stdout.decode().split("\n") if line]
    with open(expected, "rt") as exp_res_fh:
        exp_res = [line.strip().split("\t") for line in exp_res_fh.readlines() if line]
        assert res == exp_res

import os
import subprocess

import imagehash
import pytest
from PIL import Image


@pytest.mark.parametrize(
    ["bam", "bed", "expected", "config"],
    [
        # Standard case
        (
            "test/standard/HG00096_hifi.bam",
            "test/standard/region.bed",
            "test/standard/expected.bed",
            tuple(["-c", "test/standard/config.toml"]),
        ),
        # Ignore regions
        (
            "test/ignored/HG00731_hifi.bam",
            "test/ignored/region.bed",
            "test/ignored/expected.bed",
            tuple(
                [
                    "-c",
                    "test/ignored/config.toml",
                    "--ignore_regions",
                    "test/ignored/ignore.bed",
                ]
            ),
        ),
        # Percent misjoin threshold
        (
            "test/misjoin/HG00171_hifi.bam",
            "test/misjoin/region.bed",
            "test/misjoin/expected_perc.bed",
            tuple(["-c", "test/misjoin/config_perc.toml"]),
        ),
        # No reads covering region.
        (
            "test/all_gap/cov.tsv.gz",
            "test/all_gap/region.bed",
            "test/all_gap/expected.bed",
            tuple(),
        ),
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


# Check that providing no bai produces non-zero exit code.
def test_bam_idx_check():
    infile = "test/no_bai/null.bam"
    process = subprocess.run(
        ["python", "-m", "nucflag.main", "-i", infile],
        capture_output=True,
    )
    assert (
        process.returncode == 1
        and f"FileNotFoundError: {infile} must be indexed. Run 'samtools index {infile}'."
        in process.stderr.decode()
    )


@pytest.mark.parametrize(
    ["cov", "bed", "expected_dir", "output_dir", "config"],
    [
        # Overlay many beds.
        *[
            (
                "test/overlay/NA20847_rc-chr3_haplotype2-0000105:89881870-96384969.bed.gz",
                "test/overlay/region.bed",
                f"test/overlay/expected/{i}",
                f"test/overlay/output_{i}/",
                tuple(
                    [
                        "-c",
                        "test/overlay/config.toml",
                        "--overlay_regions",
                        *["test/overlay/repeatmasker.bed" for _ in range(i)],
                    ]
                ),
            )
            for i in range(1, 4)
        ],
        # Ignore a position.
        (
            "test/overlay/NA20847_rc-chr3_haplotype2-0000105:89881870-96384969.bed.gz",
            "test/overlay/region.bed",
            "test/overlay/expected/ignore/",
            "test/overlay/output_ignore/",
            tuple(
                [
                    "-c",
                    "test/overlay/config.toml",
                    "--overlay_regions",
                    "test/overlay/repeatmasker_ignore.bed",
                ]
            ),
        ),
        # Misassembly overlaps partially with normal region.
        (
            "test/overlay/AG16778_chr4_contig-0003083:3247326-8431235.bed.gz",
            "test/overlay/region_overlap_partial.bed",
            "test/overlay/expected/overlap_partial/",
            "test/overlay/overlap_partial/",
            tuple(
                [
                    "-c",
                    "test/overlay/config.toml",
                    "--ignore_regions",
                    "test/overlay/repeatmasker_overlap_partial.bed",
                ]
            ),
        ),
    ],
)
def test_correct_plot(
    cov: str, bed: str, expected_dir: str, output_dir: str, config: tuple[str]
):
    contigs = []
    with open(bed, "rt") as fh:
        for line in fh.readlines():
            name, start, stop, *_ = line.strip().split("\t")
            contigs.append(f"{name}:{start}-{stop}")

    _ = subprocess.run(
        [
            "python",
            "-m",
            "nucflag.main",
            "-i",
            cov,
            "-b",
            bed,
            "-d",
            output_dir,
            *config,
        ],
        capture_output=True,
        check=True,
    )

    for ctg in contigs:
        exp_plot_path = os.path.join(expected_dir, f"{ctg}.png")
        out_plot_path = os.path.join(output_dir, f"{ctg}.png")

        # https://stackoverflow.com/q/49595541
        # Color hash to compare features.
        assert imagehash.colorhash(Image.open(exp_plot_path)) == imagehash.colorhash(
            Image.open(out_plot_path)
        )

        # Remove outplot
        os.remove(out_plot_path)

    # Remove output dir.
    os.rmdir(output_dir)

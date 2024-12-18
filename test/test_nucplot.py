import os
import subprocess

import imagehash
import pytest
from PIL import Image

from .helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["bam", "bed", "expected_misassemblies", "expected_statuses", "config"],
    [
        # Standard case
        (
            "test/standard/HG00096_hifi.bam",
            "test/standard/region.bed",
            "test/standard/expected.bed",
            None,
            tuple(["-c", "test/standard/config.toml"]),
        ),
        # Ignore regions
        (
            "test/ignored/HG00731_hifi.bam",
            "test/ignored/region.bed",
            "test/ignored/expected.bed",
            None,
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
            None,
            tuple(["-c", "test/misjoin/config_perc.toml"]),
        ),
        # No reads covering region.
        (
            "test/all_gap/cov.tsv.gz",
            "test/all_gap/region.bed",
            "test/all_gap/expected.bed",
            None,
            tuple(),
        ),
        # Check that HETs don't affect status.
        (
            "test/status_ignore_het/HG00514_rc-chr8_haplotype1-0000015:40851933-44831382.bed.gz",
            "test/status_ignore_het/region.bed",
            "test/status_ignore_het/expected_misassemblies.bed",
            "test/status_ignore_het/expected_status.bed",
            tuple(
                [
                    "-c",
                    "test/status_ignore_het/config.toml",
                    "--ignore_regions",
                    "test/status_ignore_het/ignore.bed",
                ]
            ),
        ),
        # Check that collapses are trimmed if there is an overlapping misjoin
        (
            "test/trim_collapse_misjoin/HG02953_chr4_h1tg000006l#1-193384017:49454294-55363494:49454294-55363494.bed.gz",
            "test/trim_collapse_misjoin/region.bed",
            "test/trim_collapse_misjoin/expected_misassemblies.bed",
            "test/trim_collapse_misjoin/expected_status.bed",
            tuple(["-c", "test/trim_collapse_misjoin/config.toml"]),
        ),
        # Check that misjoins that reach 1 are detected.
        (
            "test/missing_valleys/NA21487_chr22_h1tg000028l#1-46561102:8111393-11481292.bed.gz",
            "test/missing_valleys/region.bed",
            "test/missing_valleys/expected_misassemblies.bed",
            "test/missing_valleys/expected_status.bed",
            tuple(["-c", "test/missing_valleys/config.toml"]),
        ),
    ],
)
def test_identify_misassemblies(
    bam: str,
    bed: str,
    expected_misassemblies: str,
    expected_statuses: str | None,
    config: tuple[str],
):
    expected_outputs = [("-o", expected_misassemblies)]
    if expected_statuses:
        expected_outputs.append(("-s", expected_statuses))

    run_integration_test(
        "python",
        "-m",
        "nucflag.main",
        "-i",
        bam,
        "-b",
        bed,
        *config,
        expected_output=expected_outputs,
    )


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

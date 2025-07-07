# import os
# import subprocess

# import imagehash
# import pytest
# from PIL import Image

# from .helpers.integration import run_integration_test


# # Check that providing no bai produces non-zero exit code.
# def test_bam_idx_check():
#     infile = "test/no_bai/null.bam"
#     process = subprocess.run(
#         ["python", "-m", "nucflag.main", "-i", infile],
#         capture_output=True,
#     )
#     assert (
#         process.returncode == 1
#         and f"FileNotFoundError: {infile} must be indexed. Run 'samtools index {infile}'."
#         in process.stderr.decode()
#     )


# @pytest.mark.parametrize(
#     ["cov", "bed", "expected_dir", "output_dir", "config"],
#     [
#         # Overlay many beds.
#         *[
#             (
#                 "test/overlay/NA20847_rc-chr3_haplotype2-0000105:89881870-96384969.bed.gz",
#                 "test/overlay/region.bed",
#                 f"test/overlay/expected/{i}",
#                 f"test/overlay/output_{i}/",
#                 tuple(
#                     [
#                         "-c",
#                         "test/overlay/config.toml",
#                         "--overlay_regions",
#                         *["test/overlay/repeatmasker.bed" for _ in range(i)],
#                     ]
#                 ),
#             )
#             for i in range(1, 4)
#         ],
#         # Ignore a position.
#         (
#             "test/overlay/NA20847_rc-chr3_haplotype2-0000105:89881870-96384969.bed.gz",
#             "test/overlay/region.bed",
#             "test/overlay/expected/ignore/",
#             "test/overlay/output_ignore/",
#             tuple(
#                 [
#                     "-c",
#                     "test/overlay/config.toml",
#                     "--overlay_regions",
#                     "test/overlay/repeatmasker_ignore.bed",
#                 ]
#             ),
#         ),
#         # Misassembly overlaps partially with normal region.
#         (
#             "test/overlay/AG16778_chr4_contig-0003083:3247326-8431235.bed.gz",
#             "test/overlay/region_overlap_partial.bed",
#             "test/overlay/expected/overlap_partial/",
#             "test/overlay/overlap_partial/",
#             tuple(
#                 [
#                     "-c",
#                     "test/overlay/config_partial.toml",
#                     "--ignore_regions",
#                     "test/overlay/repeatmasker_overlap_partial.bed",
#                 ]
#             ),
#         ),
#     ],
# )
# def test_correct_plot(
#     cov: str, bed: str, expected_dir: str, output_dir: str, config: tuple[str]
# ):
#     contigs = []
#     with open(bed, "rt") as fh:
#         for line in fh.readlines():
#             name, start, stop, *_ = line.strip().split("\t")
#             contigs.append(f"{name}:{start}-{stop}")

#     _ = subprocess.run(
#         [
#             "python",
#             "-m",
#             "nucflag.main",
#             "-i",
#             cov,
#             "-b",
#             bed,
#             "-d",
#             output_dir,
#             *config,
#         ],
#         capture_output=True,
#         check=True,
#     )

#     for ctg in contigs:
#         exp_plot_path = os.path.join(expected_dir, f"{ctg}.png")
#         out_plot_path = os.path.join(output_dir, f"{ctg}.png")

#         # https://stackoverflow.com/q/49595541
#         # Color hash to compare features.
#         assert imagehash.colorhash(Image.open(exp_plot_path)) == imagehash.colorhash(
#             Image.open(out_plot_path)
#         )

#         # Remove outplot
#         os.remove(out_plot_path)

#     # Remove output dir.
#     os.rmdir(output_dir)

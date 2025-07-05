import os


def check_indexed(infile: str) -> None:
    if infile.endswith(".bam") and not (
        os.path.exists(f"{infile}.bai") or os.path.exists(f"{infile}.csi")
    ):
        raise FileNotFoundError(
            f"{infile} must be indexed. Run 'samtools index {infile}'."
        )
    elif infile.endswith(".cram") and not os.path.exists(f"{infile}.crai"):
        raise FileNotFoundError(
            f"{infile} must be indexed. Run 'samtools index {infile}'."
        )

import os


def check_bam_indexed(infile: str) -> None:
    if infile.endswith(".bam") and not os.path.exists(f"{infile}.bai"):
        raise FileNotFoundError(
            f"{infile} must be indexed. Run 'samtools index {infile}'."
        )

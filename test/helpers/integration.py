import sys
import gzip
import subprocess
import itertools
import tempfile
from typing import Iterable

Output_Lines = list[str]
Outputs = Output_Lines | str
Outputs_To_Check = list[tuple[Outputs, str]]


def check_output(outputs: Outputs_To_Check, overwrite_output: bool) -> None:
    """
    Check that outputs match expected output line-by-line.
    """
    for in_output, exp_output in outputs:
        if isinstance(in_output, list):
            sorted_in_res = sorted(in_output)
        else:
            with open(in_output, "rt") as fh:
                sorted_in_res = sorted(line.strip() for line in fh.readlines())

        if overwrite_output:
            is_gzip_file = exp_output.endswith(".gz")
            with gzip.open(exp_output, "wb") if is_gzip_file else open(
                exp_output, "wt"
            ) as fh:
                for line in sorted_in_res:
                    fh.write(f"{line}\n".encode() if is_gzip_file else f"{line}\n")
            continue

        with gzip.open(exp_output, "rt") if exp_output.endswith(".gz") else open(
            exp_output, "rt"
        ) as exp_res_fh:
            sorted_exp_res = sorted(
                line.strip() for line in exp_res_fh.readlines() if line
            )
            for i, (res, exp) in enumerate(
                zip(sorted_in_res, sorted_exp_res, strict=True)
            ):
                assert res == exp, f"Not equal ({res} != {exp}) on line {i}"


def run_integration_test(
    *cmd: str, expected_output: str | Iterable[tuple[str, str]], overwrite_output: bool
) -> None:
    """
    Run integration test and check/cleans up outputs.

    # Args
    * `cmd`
        * Command to run.
        * If `expected_output` is an iterable, the outputs are appended to the end of the cmd.
    * `expected_output`
        * Either a single output or multiple outputs.
        * If the output is an iterator, expects a 2-element tuple:
            1. output command flag
            2. expected output.
    * `overwrite_output`
        * For each output, write result to output file. Use to regenerate expected files.
    """
    if isinstance(expected_output, str):
        process = subprocess.run(
            [*cmd],
            capture_output=True,
            check=True,
        )
        res = sorted(
            line.strip() for line in process.stdout.decode().split("\n") if line
        )
        outputs: Outputs_To_Check = [(res, expected_output)]
        check_output(outputs, overwrite_output)
    else:
        flags, expected = zip(*expected_output)
        outfiles = [tempfile.NamedTemporaryFile() for _ in flags]
        outfile_names = [file.name for file in outfiles]
        new_cmd = [*cmd, *list(itertools.chain(*zip(flags, outfile_names)))]
        print(" ".join(new_cmd), file=sys.stderr)
        subprocess.run(
            [*new_cmd],
            check=True,
        )
        outputs = [out for out in zip(outfile_names, expected)]
        check_output(outputs, overwrite_output)
        for file in outfiles:
            file.close()

import sys
import subprocess
import itertools
import tempfile
import polars as pl
from typing import Iterable
from polars.testing import assert_frame_equal


Outputs_To_Check = list[tuple[pl.DataFrame | str, str]]


def check_output(outputs: Outputs_To_Check, overwrite_output: bool) -> None:
    """
    Check that outputs match expected output line-by-line.
    """
    for in_output, exp_output in outputs:
        if isinstance(in_output, pl.DataFrame):
            df_in_res = in_output
        else:
            df_in_res = pl.read_csv(in_output, separator="\t")

        if overwrite_output:
            df_in_res.write_csv(exp_output, separator="\t", include_header=True)
            continue

        df_exp_res = pl.read_csv(exp_output, separator="\t")
        assert_frame_equal(df_in_res, df_exp_res, check_row_order=False)


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
        res = pl.read_csv(process.stdout, separator="\t")
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

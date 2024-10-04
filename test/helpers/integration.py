import sys
import subprocess
import itertools
import tempfile
from typing import Iterable

Output_Lines = list[str]
Outputs = Output_Lines | str
Outputs_To_Check = list[tuple[Outputs, str]]


def check_output(outputs: Outputs_To_Check) -> None:
    """
    Check that outputs match expected output line-by-line.
    """
    for in_output, exp_output in outputs:
        if isinstance(in_output, list):
            sorted_in_res = sorted(in_output)
        else:
            with open(in_output, "rt") as fh:
                sorted_in_res = sorted(line.strip() for line in fh.readlines())

        with open(exp_output, "rt") as exp_res_fh:
            sorted_exp_res = sorted(
                line.strip() for line in exp_res_fh.readlines() if line
            )
            assert sorted_in_res == sorted_exp_res


def run_integration_test(
    *cmd: str, expected_output: str | Iterable[tuple[str, str]]
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
        check_output(outputs)
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
        check_output(outputs)
        for file in outfiles:
            file.close()

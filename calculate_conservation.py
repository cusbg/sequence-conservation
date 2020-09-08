#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import tempfile
import shutil
import typing
import logging
import time
import argparse
import json

MARK_SEQUENCE_PREFIX = "> query_sequence |"

PSIBLAST_CMD = "/opt/conservation/ncbi-blast-2.9.0+/bin/psiblast"
BLASTDBCMD_CMD = "/opt/conservation/ncbi-blast-2.9.0+/bin/blastdbcmd"
CDHIT_CMD = "/opt/conservation/cd-hit-v4.8.1-2019-0228/cd-hit"
MUSCLE_CMD = "/opt/conservation/muscle3.8.31_i86linux64"
JENSE_SHANNON_DIVERGANCE_DIR = "/opt/conservation/conservation_code/"

MIN_SEQUENCE_COUNT = 50

time_end_before = None


def _read_arguments() -> typing.Dict[str, str]:
    parser = argparse.ArgumentParser(
        description="Compute conservation scores for given sequences.")
    parser.add_argument(
        "--input", required=True,
        help="Input FASTA file.")
    parser.add_argument(
        "--output", required=True,
        help="Output JSON-lines file.")
    parser.add_argument(
        "--output-raw", default=None,
        help="Output file with raw data with '%s' "
             "to be substituted with chain index.")
    parser.add_argument(
        "--time-limit", default=None, type=int,
        help="Soft time limit in seconds for running external tools.")
    return vars(parser.parse_args())


def main(arguments):
    _init_logging()

    time_start = time.time()
    _set_time_out(arguments)

    working_root_dir = tempfile.mkdtemp("", "conservation-")
    logging.info("Processing file: %s", arguments["input"])

    if _output_exists(arguments):
        logging.info("Output file already exists.")
        return

    try:
        sequences = _iterate_fasta_file(
            arguments["input"], lambda line: line.rstrip())
        process_sequences(sequences, working_root_dir, arguments)
    finally:
        logging.info("Execution time: %1.2f s", time.time() - time_start)
        shutil.rmtree(working_root_dir, ignore_errors=True)


def _output_exists(arguments) -> bool:
    if not os.path.exists(arguments["output"]):
        return False
    if arguments["output-raw"] is not None and \
            not os.path.exists(arguments["output-raw"]):
        return False
    return True


def _init_logging() -> None:
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S")


def _set_time_out(arguments):
    if arguments["time_limit"] is None:
        return
    global time_end_before
    time_end_before = time.time() + arguments["time_limit"]


def _iterate_fasta_file(input_file: str, on_line=lambda line: line) \
        -> typing.Iterable[typing.Tuple[str, str]]:
    header = None
    sequence = ""
    with open(input_file) as in_stream:
        for line in in_stream:
            line = on_line(line)
            if line.startswith(">"):
                if header is None:
                    header = line
                else:
                    yield header, sequence
                    header = line
                    sequence = ""
            else:
                sequence += line
    if header is not None:
        yield header, sequence


def process_sequences(
        sequences: typing.Iterable[typing.Tuple[str, str]],
        working_root_dir: str, arguments) -> None:
    global time_end_before
    # We write to temp file and then move only once done.
    temp_output = arguments["output"] + ".tmp"
    with open(temp_output, "w", encoding="utf-8") as out_steam:
        for index, (header, sequence) in enumerate(sequences):
            if time_end_before is not None and time_end_before < time.time():
                raise TimeoutError()
            working_dir = \
                os.path.join(working_root_dir, str(index).zfill(6))
            os.makedirs(working_dir)
            result = compute_conservation(sequence, header, working_dir)
            if result["file"] is None:
                _on_fail_to_compute(header, out_steam)
                continue
            # Save result into JSON.
            json.dump({
                "header": header,
                "database": result["database"],
                "conservation": _load_conservation_result_file(result["file"])
            }, out_steam)
            out_steam.write("\n")
            # Save RAW file.
            if arguments.get("output-raw", None) is not None:
                os.rename(result["file"], arguments["output-raw"].format(index))

    os.rename(temp_output, arguments["output"])


def compute_conservation(
        sequence: str, header: str, working_dir: str) -> typing.Dict:
    pdb_file = os.path.join(working_dir, "input-sequence.fasta")
    _save_sequence_to_pdb(MARK_SEQUENCE_PREFIX + header, sequence, pdb_file)

    blast_output_file = os.path.join(working_dir, "blast-output")
    _blast_sequence(pdb_file, "swissprot", blast_output_file, working_dir)
    used_database = "swissprot"
    if not _enough_blast_results(blast_output_file, MIN_SEQUENCE_COUNT):
        _blast_sequence(pdb_file, "uniref90", blast_output_file, working_dir)
        used_database = "uniref90"
        if not _enough_blast_results(blast_output_file, MIN_SEQUENCE_COUNT):
            logging.error("Not enough sequences found!")
            return {
                "file": None,
                "database": None
            }

    muscle_output_file = os.path.join(working_dir, "muscle-output")
    _execute_muscle(
        pdb_file, blast_output_file, muscle_output_file, working_dir)

    conservation_input_file = os.path.join(working_dir, "conservation-input")
    _order_muscle_result(muscle_output_file, conservation_input_file)

    conservation_output_file = os.path.join(working_dir, "conservation-output")
    return {
        "file": _execute_jensen_shannon_divergence(
            conservation_input_file, conservation_output_file),
        "database": used_database
    }


def _blast_sequence(
        pdb_file: str, database: str, output_file: str, working_dir: str) \
        -> None:
    # Find similar sequences.
    logging.info("Running psiblast on '%s' database ...", database)
    psiblast_file = os.path.join(working_dir, "psiblast")
    _execute_psiblast(pdb_file, psiblast_file, database)
    # Filter results.
    logging.info("Filtering files ...")
    psiblast_file_filtered = os.path.join(working_dir, "psiblast-filtered")
    _filter_psiblast_file(psiblast_file, psiblast_file_filtered)
    # Get sequences for results from previous step.
    logging.info("Running blastdbcmd ...")
    sequences_file = os.path.join(working_dir, "blastdb-output")
    _execute_blastdbcmd(psiblast_file_filtered, sequences_file, database)
    # Cluster and select representatives.
    logging.info("Running cd-hit ...")
    cdhit_log_file = os.path.join(working_dir, "cd-hit.log")
    _execute_cdhit(sequences_file, output_file, cdhit_log_file)


def _save_sequence_to_pdb(header: str, sequence: str, output_file: str) -> None:
    line_width = 80
    lines = [
        sequence[index:index + line_width]
        for index in range(0, len(sequence), line_width)
    ]
    with open(output_file, "w") as out_stream:
        out_stream.write(header)
        out_stream.write("\n")
        out_stream.write("\n".join(lines))
        out_stream.write("\n")


def _execute_psiblast(pdb_file: str, output_file: str, database: str) -> None:
    output_format = "6 sallseqid qcovs pident"
    cmd = "{} < {} -db {} -outfmt '{}' -evalue 1e-5 > {}".format(
        PSIBLAST_CMD, pdb_file, database, output_format, output_file)
    logging.debug("Executing BLAST ...")
    _execute(cmd)


def _execute(command: str):
    global time_end_before
    if time_end_before is None:
        logging.debug("Executing command:\n%s", command)
        subprocess.run(command, shell=True, env=os.environ.copy())
        return
    if time_end_before < time.time():
        raise TimeoutError()
    timeout = time_end_before - time.time()
    logging.debug("Executing with timeout: %s command:\n%s", timeout, command)
    subprocess.run(command, shell=True, env=os.environ.copy(), timeout=timeout)


def _filter_psiblast_file(input_file: str, output_file: str):
    results_count = 0
    with open(input_file) as in_stream:
        with open(output_file, "w") as out_stream:
            for line in in_stream:
                identifier, coverage, identity = line.rstrip().split("\t")
                if _filter_condition(float(coverage), float(identity)):
                    out_stream.write(identifier)
                    out_stream.write("\n")
                    results_count += 1
    return results_count


def _filter_condition(coverage: float, identity: float) -> bool:
    return coverage >= 80 and 30 <= identity <= 95


def _execute_blastdbcmd(
        psiblast_output_file: str, sequence_file: str, database: str) -> None:
    cmd = "{} -db {} -entry_batch {} > {}".format(
        BLASTDBCMD_CMD, database, psiblast_output_file,
        sequence_file)
    logging.debug("Executing BLAST ...")
    _execute(cmd)


def _execute_cdhit(input_file: str, output_file: str, log_file: str) -> None:
    cmd = "{} -i {} -o {} > {}".format(
        CDHIT_CMD, input_file, output_file, log_file)
    logging.debug("Executing CD-HIT ..")
    _execute(cmd)


def _enough_blast_results(fasta_file: str, min_count: int) -> bool:
    counter = 0
    for _, _ in _iterate_fasta_file(fasta_file):
        counter += 1
    logging.info("Number of sequences in %s is %s", fasta_file, counter)
    return counter > min_count


def _execute_muscle(
        pdb_file: str, sequence_file: str, output_file: str, working_dir: str) \
        -> None:
    muscle_input = os.path.join(working_dir, "muscle-input")
    _merge_files([sequence_file, pdb_file], muscle_input)
    cmd = "cat {} | {} -quiet > {}".format(
        muscle_input, MUSCLE_CMD, output_file)
    logging.info("Executing muscle ...")
    _execute(cmd)


def _merge_files(input_files: typing.List[str], output_file: str) -> None:
    with open(output_file, "w") as out_stream:
        for input_file in input_files:
            with open(input_file) as in_stream:
                for line in in_stream:
                    out_stream.write(line)


def _order_muscle_result(input_file: str, output_file: str) -> None:
    """Put the marked sequence at the top of the file and change it's header."""
    logging.info("Ordering muscle results ...")
    first_header = None
    first_sequence = None
    for header, sequence in _iterate_fasta_file(input_file):
        if header.startswith(MARK_SEQUENCE_PREFIX):
            # We can remove the prefix here
            first_header = header[len(MARK_SEQUENCE_PREFIX):]
            first_sequence = sequence
            break
    if first_header is None:
        raise Exception("Missing header '" + MARK_SEQUENCE_PREFIX +
                        "' in " + input_file)
    with open(output_file, "w") as out_stream:
        out_stream.write(first_header)
        out_stream.write(first_sequence)
        for header, sequence in _iterate_fasta_file(input_file):
            if header.startswith(MARK_SEQUENCE_PREFIX):
                continue
            out_stream.write(header)
            out_stream.write(sequence)


def _execute_jensen_shannon_divergence(
        input_file: str, output_file: str) -> str:
    """Input sequence must be on the first position."""
    cmd = "cd {} && python2 score_conservation.py {} > {}".format(
        JENSE_SHANNON_DIVERGANCE_DIR,
        os.path.abspath(input_file),
        os.path.abspath(output_file))
    logging.info("Executing Jense Shannon Divergence script ...")
    logging.debug("Executing command:\n%s", cmd)
    _execute(cmd)
    return output_file


def _on_fail_to_compute(header: str, out_steam) -> None:
    json.dump({
        "header": header,
        "error": "Not enough BLAST results."
    }, out_steam)
    out_steam.write("\n")


def _load_conservation_result_file(conservation_output: str) \
        -> typing.List[float]:
    """
    Return conservation for first sequence in input.
    We know that we put our sequence on the first place, so we can use
    this information to retrieve it's conservation.
    We use lines where the third column contains AA, i.e. does not start
    with '-'.
    """
    result = []
    with open(conservation_output) as in_stream:
        # Read header.
        for line in in_stream:
            if not line.startswith("#"):
                break
        # Read data.
        for line in in_stream:
            tokens = [token for token in line.rstrip().split("\t")]
            if tokens[2].startswith("-"):
                continue
            result.append(float(tokens[1]))
    return result


if __name__ == "__main__":
    main(_read_arguments())

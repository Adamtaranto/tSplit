"""
Utility functions for tSplit.

This module provides helper functions for file operations, dependency checking,
sequence ID handling, and other utility tasks used across the TE-splitter application.
"""

import logging
import os
import re
import shutil
import sys
from typing import Any, Generator, List, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def check_tools(
    required_tools: Optional[List[str]] = None,
    optional_tools: Optional[List[str]] = None,
) -> None:
    """
    Check if required and optional tools are available on the system's PATH.

    Parameters
    ----------
    required_tools : list of str, optional
        List of required tool names.
    optional_tools : list of str, optional
        List of optional tool names.

    Raises
    ------
    RuntimeError
        If any required tool is not found.

    Notes
    -----
    Tools in the required_tools list will cause the program to exit if not found.
    Tools in the optional_tools list will only produce a warning if not found.
    """
    # Initialize empty lists if None is provided
    if required_tools is None:
        required_tools = []
    if optional_tools is None:
        optional_tools = []

    missing_required_tools = []

    def print_message(tool: str, path: Optional[str], color_code: str) -> None:
        """
        Print a message to stderr with the tool name and path in the specified color.

        Parameters
        ----------
        tool : str
            The name of the tool.
        path : str or None
            The path to the tool, or None if not found.
        color_code : str
            The ANSI color code for the message.
        """
        tool_padded = tool.ljust(15)
        if path:
            message = f'{color_code}{tool_padded}\t{path}\033[0m'
        else:
            message = f'{color_code}{tool_padded}\tNOT FOUND\033[0m'
        print(message, file=sys.stderr)

    # Check required tools
    logging.info('Checking for dependencies:')
    for tool in required_tools:
        path = shutil.which(tool)  # Get the full path of the executable
        if path:
            print_message(tool, path, '\033[92m')  # Green for found tools
        else:
            print_message(tool, None, '\033[91m')  # Red for missing required tools
            missing_required_tools.append(tool)

    # Check optional tools
    for tool in optional_tools:
        path = shutil.which(tool)
        if path:
            print_message(tool, path, '\033[92m')  # Green for found tools
        else:
            print_message(tool, None, '\033[93m')  # Yellow for missing optional tools

    # Raise error if any required tool is missing
    if missing_required_tools:
        error_message = 'ERROR: Some required tools could not be found: ' + ', '.join(
            missing_required_tools
        )
        logging.error(error_message)
        raise RuntimeError(
            'Missing required tools: ' + ', '.join(missing_required_tools)
        )


def tSplitchecks(args: Any) -> str:
    """
    Perform housekeeping tasks for TE-splitter.

    Creates output directories, sets file naming conventions,
    and validates input files.

    Parameters
    ----------
    args : argparse.Namespace
        Argument parser object containing input parameters.

    Returns
    -------
    str
        Full path to the output file.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    """
    # Check if the input file exists
    if not os.path.isfile(args.infile):
        logging.error('Input sequence file does not exist. Quitting.')
        raise FileNotFoundError(f"Input sequence file '{args.infile}' does not exist.")
    logging.info(f'Input file found: {args.infile}')

    # Determine the output directory - create if it doesn't exist
    if args.outdir:
        absOutDir = os.path.abspath(args.outdir)
        if not os.path.isdir(absOutDir):
            os.makedirs(absOutDir)
            logging.info(f'Creating output directory: {absOutDir}')
        outDir = absOutDir
    else:
        outDir = os.getcwd()  # Use current working directory if none specified
    logging.debug(f'Set output directory: {outDir}')

    # Set the prefix for output files - use input filename if none provided
    if not args.prefix:
        prefix = os.path.splitext(os.path.basename(args.infile))[0]
    else:
        prefix = args.prefix
    logging.debug(f'Set prefix: {prefix}')

    # Create the output file path by combining directory and filename
    outfile = prefix + '_tsplit_output.fasta'
    outpath = os.path.join(outDir, outfile)
    logging.debug(f'Set outfile target: {outpath}')

    # Return the full path to the output file
    return outpath


def cleanID(s: str) -> str:
    """
    Remove non-alphanumeric characters from string and replace whitespace.

    Parameters
    ----------
    s : str
        Input string to be cleaned.

    Returns
    -------
    str
        Cleaned string with special characters removed and
        whitespace replaced with underscores.

    Notes
    -----
    This function is useful for ensuring sequence IDs are compatible
    with various bioinformatics tools.
    """
    # Remove any character that isn't alphanumeric or whitespace
    s = re.sub(r'[^\w\s]', '', s)
    # Replace any whitespace sequence with a single underscore
    s = re.sub(r'\s+', '_', s)
    return s


def segWrite(
    outfile: str, segs: Optional[Generator[SeqRecord, None, None]] = None
) -> None:
    """
    Write sequence records to a FASTA file.

    Takes a generator object yielding SeqRecord objects and
    writes each to the specified output file in FASTA format.
    If no sequences are written, the output file is removed.

    Parameters
    ----------
    outfile : str
        Path to the output file.
    segs : generator of Bio.SeqRecord.SeqRecord, optional
        Generator yielding SeqRecord objects to write.

    Returns
    -------
    None
        This function writes data to a file and doesn't return any value.

    Notes
    -----
    If no sequences are written (empty generator), the output file
    will be automatically removed to avoid empty files.
    """
    seqcount = 0
    if segs:
        # Open the output file for writing
        with open(outfile, 'w') as handle:
            for seq in segs:
                seqcount += 1
                SeqIO.write(seq, handle, 'fasta')

        # Clean up empty files to avoid clutter
        if seqcount == 0:
            os.remove(outfile)


def write_paf(alignment_data: List[dict], outfile: str) -> None:
    """
    Write alignment data in PAF (Pairwise mApping Format).

    PAF is a text format for representing alignment between two sequences.
    Each line represents a single alignment with 12 mandatory fields.

    Parameters
    ----------
    alignment_data : list of dict
        List of dictionaries containing alignment information.
        Each dict should have keys: qry_name, qry_length, qry_start, qry_end,
        strand, ref_name, ref_length, ref_start, ref_end, num_matches,
        aln_block_length, mapping_quality.
    outfile : str
        Path to the output PAF file.

    Returns
    -------
    None
        This function writes data to a file and doesn't return any value.

    Notes
    -----
    PAF format specification:
    1. Query sequence name
    2. Query sequence length
    3. Query start (0-based, inclusive)
    4. Query end (0-based, exclusive)
    5. Relative strand: "+" or "-"
    6. Target sequence name
    7. Target sequence length
    8. Target start on original strand (0-based, inclusive)
    9. Target end on original strand (0-based, exclusive)
    10. Number of matching bases
    11. Alignment block length
    12. Mapping quality (0-255; 255 for missing)

    Positions use half-open intervals [start, end) where start is inclusive
    and end is exclusive (0-based coordinates).
    """
    logging.info(f'Writing PAF output to: {outfile}')
    with open(outfile, 'w') as f:
        for aln in alignment_data:
            # Write PAF line with 12 mandatory fields
            f.write(
                f'{aln["qry_name"]}\t{aln["qry_length"]}\t{aln["qry_start"]}\t{aln["qry_end"]}\t'
                f'{aln["strand"]}\t{aln["ref_name"]}\t{aln["ref_length"]}\t{aln["ref_start"]}\t'
                f'{aln["ref_end"]}\t{aln["num_matches"]}\t{aln["aln_block_length"]}\t'
                f'{aln["mapping_quality"]}\n'
            )
    logging.info(f'Wrote {len(alignment_data)} alignments to PAF file.')


def write_gff3(feature_data: List[dict], outfile: str, source: str = 'tSplit') -> None:
    """
    Write terminal repeat features in GFF3 format.

    GFF3 (General Feature Format version 3) is a standard format for
    representing genomic features.

    Parameters
    ----------
    feature_data : list of dict
        List of dictionaries containing feature information.
        Each dict should have keys: seqid, type, start, end, strand, attributes.
        Optional keys: score, phase.
    outfile : str
        Path to the output GFF3 file.
    source : str, optional
        Source of the annotations (default: 'tSplit').

    Returns
    -------
    None
        This function writes data to a file and doesn't return any value.

    Notes
    -----
    GFF3 format has 9 tab-separated columns:
    1. seqid - sequence ID
    2. source - source of the annotation (e.g., 'tSplit')
    3. type - feature type (e.g., 'TIR', 'LTR')
    4. start - 1-based start position
    5. end - 1-based end position
    6. score - score (or '.' for none)
    7. strand - '+', '-', or '.' for unknown
    8. phase - CDS phase (or '.' for non-CDS features)
    9. attributes - semicolon-separated attribute=value pairs
    """
    logging.info(f'Writing GFF3 output to: {outfile}')
    with open(outfile, 'w') as f:
        # Write GFF3 header
        f.write('##gff-version 3\n')

        # Write each feature
        for feat in feature_data:
            score = feat.get('score', '.')
            phase = feat.get('phase', '.')
            # Write GFF3 line with 9 columns
            f.write(
                f'{feat["seqid"]}\t{source}\t{feat["type"]}\t{feat["start"]}\t{feat["end"]}\t'
                f'{score}\t{feat["strand"]}\t{phase}\t{feat["attributes"]}\n'
            )
    logging.info(f'Wrote {len(feature_data)} features to GFF3 file.')

"""
TIR identification and extraction command module.

This module implements the command-line functionality for identifying and
extracting Terminal Inverted Repeats (TIRs) from DNA transposons. It handles
the initialization of the TIR analysis pipeline, including setting up logging,
checking dependencies, and orchestrating the extraction and output of TIRs.

The module works with the parseAlign module to perform self-alignment of
sequences and identify terminal repeats that occur in opposite orientations
(characteristic of TIRs in DNA transposons). It can also optionally construct
synthetic Miniature Inverted-repeat Transposable Elements (MITEs).
"""

from argparse import Namespace
from typing import Optional

from tsplit.logs import init_logging
from tsplit.parseAlign import getTIRs
from tsplit.utils import check_tools, segWrite, tSplitchecks


def main(args: Optional[Namespace] = None) -> None:
    """
    Execute the TIR identification and extraction workflow.

    Orchestrates the process of identifying Terminal Inverted Repeats in
    DNA transposons by setting up logging, checking tool dependencies,
    preparing output paths, running the TIR analysis, and writing results.

    Parameters
    ----------
    args : argparse.Namespace, optional
        Command line arguments containing configuration for TIR analysis.
        If None, the function will use previously parsed arguments.

    Returns
    -------
    None
        Results are written to output files specified in args.

    Notes
    -----
    This function requires external tools for sequence alignment:
    - Required: delta-filter, nucmer, show-coords
    - Optional or required based on method: blastn

    The function uses the getTIRs function from parseAlign module to identify
    TIRs and extract segments based on the specified reporting mode. If the
    'makemites' option is enabled, synthetic MITEs will be constructed from
    the identified TIRs.
    """
    # Set up logging with specified verbosity level
    init_logging(loglevel=args.loglevel)

    # Check for required external programs
    # For nucmer alignment, we need the MUMmer suite tools
    required_tools = ['delta-filter', 'nucmer', 'show-coords']
    optional_tools = ['blastn']

    # If using BLAST for alignment instead of nucmer, adjust requirements
    if args.method == 'blastn':
        required_tools.append('blastn')
        optional_tools = []  # blastn is now required, not optional

    # Verify all required tools are available on the system
    check_tools(required_tools=required_tools, optional_tools=optional_tools)

    # Create output directories and validate input files
    outpath = tSplitchecks(args)

    # Search for inverted terminal repeats
    # Optionally construct synthetic MITE from TIRs if detected
    segments = getTIRs(
        args.infile,
        flankdist=args.maxdist,  # Maximum distance from element boundaries
        minterm=args.minterm,  # Minimum length of terminal repeats
        minseed=args.minseed,  # Minimum seed length for alignment
        minid=args.minid,  # Minimum percent identity between repeats
        diagfactor=args.diagfactor,  # Diagonal factor for clustering matches
        mites=args.makemites,  # Whether to construct synthetic MITEs
        report=args.splitmode,  # Mode for reporting results (split, internal, etc.)
        alignTool=args.method,  # Alignment tool (blastn or nucmer)
        temp=args.outdir,  # Directory for temporary files
        keeptemp=args.keeptemp,  # Whether to keep temporary files
    )

    # Write the identified segments to output file
    segWrite(outpath, segs=segments)

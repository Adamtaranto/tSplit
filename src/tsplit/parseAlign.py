"""
Transposable element alignment and feature extraction module.

This module provides functionality for:
1. Self-alignment of DNA sequences to identify terminal repeats
2. Extraction of Terminal Inverted Repeats (TIRs) from DNA transposons
3. Extraction of Long Terminal Repeats (LTRs) from retrotransposons
4. Construction of synthetic Miniature Inverted-repeat Transposable Elements (MITEs)
5. Filtering alignment results based on specific criteria for TEs
"""

import logging
import os
import shutil
import tempfile
from typing import Generator, List, Optional

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord
from pymummer import coords_file, nucmer

from tsplit.utils import cleanID
from tsplit.wrapping import makeBlast, run_cmd


def getTIRs(
    fasta_file: str,
    flankdist: int = 10,
    minid: float = 80,
    minterm: int = 10,
    minseed: int = 5,
    diagfactor: float = 0.3,
    mites: bool = False,
    report: str = 'split',
    temp: Optional[str] = None,
    keeptemp: bool = False,
    alignTool: str = 'nucmer',
    verbose: bool = True,
    both: bool = False,
) -> Generator[SeqRecord, None, None]:
    """
    Align elements to self and attempt to identify TIRs.

    Processes sequences to identify Terminal Inverted Repeats (TIRs)
    by performing self-alignment. Can optionally construct synthetic
    Miniature Inverted-repeat Transposable Elements (MITEs).

    Parameters
    ----------
    fasta_file : str
        Path to the multifasta file containing sequence records.
    flankdist : int, optional
        Maximum distance from element start for TIR candidates, by default 10.
    minid : float, optional
        Minimum identity between terminal repeat pairs, by default 80.
    minterm : int, optional
        Minimum length for a terminal repeat to be considered, by default 10.
    minseed : int, optional
        Minimum seed length for nucmer, by default 5.
    diagfactor : float, optional
        Diagonal factor for nucmer, by default 0.3.
    mites : bool, optional
        Whether to attempt to construct synthetic MITEs, by default False.
    report : str, optional
        Reporting mode for TIRs ('split', 'external', 'internal', 'all'), by default 'split'.
    temp : str, optional
        Path to the temporary directory, by default None.
    keeptemp : bool, optional
        Whether to keep the temporary directory after processing, by default False.
    alignTool : str, optional
        Alignment tool to use ('nucmer' or 'blastn'), by default 'nucmer'.
    verbose : bool, optional
        Whether to print verbose output, by default True.
    both : bool, optional
        Whether to report both left and right terminal repeats, by default False.

    Yields
    ------
    SeqRecord
        Segments of the sequence based on the reporting mode:
        - 'split': TIRs and internal regions separately
        - 'external': Only TIRs
        - 'internal': Only internal regions
        - 'all': Original sequences plus all segments

    Notes
    -----
    When mites=True, the function will also yield synthetic MITEs constructed by
    joining the identified TIRs.

    When both=True and report is in ['split', 'external', 'all'], both left and right
    TIRs will be yielded with suffixes '_L_TIR' and '_R_TIR'. The right TIR will be
    reverse complemented so it can be aligned with the left TIR.
    """
    # Set temp directory to cwd if none is provided
    if not temp:
        temp = os.getcwd()

    # Create a unique temporary directory
    tempDir = tempfile.mkdtemp(prefix=f'tsplit_TIR_temp_{alignTool}_', dir=temp)
    logging.debug(f'Temporary directory created: {tempDir}')

    seen_ids = set()

    found_TIRs_count = 0

    try:
        # Iterate over each record in the fasta file
        for rec in SeqIO.parse(fasta_file, 'fasta'):
            # Log the record name and length
            logging.info(
                f'Processing record {len(seen_ids) + 1}: Name: {rec.id}, Length: {len(rec)}bp'
            )

            # Check for duplicate IDs
            if rec.id in seen_ids:
                logging.error(f'Duplicate record ID found: {rec.id}')
                raise ValueError(f'Duplicate record ID found: {rec.id}')

            # Add the record ID to the set of seen IDs
            seen_ids.add(rec.id)

            # Create temp paths for single element fasta and alignment coords
            tempFasta = os.path.join(tempDir, cleanID(rec.id) + '.fasta')
            tempCoords = os.path.join(
                tempDir, cleanID(rec.id) + '_' + alignTool + '.coords'
            )

            # Write current element to single fasta
            with open(tempFasta, 'w') as f:
                SeqIO.write(rec, f, 'fasta')

            # Align to self with nucmer
            if alignTool == 'nucmer':
                # Compose Nucmer script for current element vs self
                runner = nucmer.Runner(
                    tempFasta,
                    tempFasta,
                    tempCoords,
                    min_id=minid,
                    min_length=minseed,
                    diagfactor=diagfactor,
                    mincluster=minterm,
                    breaklen=200,
                    maxmatch=True,
                    simplify=False,
                )
                # Execute nucmer
                runner.run()
            elif alignTool == 'blastn':
                # Alternatively, use blastn as search tool and write nucmer.coords-like output.
                cmd = makeBlast(seq=tempFasta, outfile=tempCoords, pid=minid)
                run_cmd(cmd, verbose=verbose, workingDir=tempDir)

            alignments = filterCoordsFileTIR(
                tempCoords, rec, minterm=minterm, flankdist=flankdist
            )

            # If alignments exist after filtering, report features using alignment
            # pair with largest internal segment i.e. first element in sorted list.
            if alignments:
                # Increment the count of elements with TIRs
                found_TIRs_count += 1

                logging.info(
                    f'Found {len(alignments)} alignments for candidate element: {rec.id}'
                )

                [print(x) for x in alignments]

                logging.info('Selecting candidate TIRs with largest internal segment:')

                print(alignments[0])

                logging.debug(
                    f'Identified coords in top candidate: \nref_start: {alignments[0].ref_start} \nref_end: {alignments[0].ref_end} \nqry_start: {alignments[0].qry_start} \nqry_end: {alignments[0].qry_end} \nref_coords: {alignments[0].ref_coords()} \nqry_coords: {alignments[0].qry_coords()} \nref_length: {alignments[0].ref_length} \nqry_length: {alignments[0].qry_length} \nframe: {alignments[0].frame} \npercent_identity: {alignments[0].percent_identity} \nhit_length_ref: {alignments[0].hit_length_ref} \nhit_length_qry: {alignments[0].hit_length_qry} \non_same_strand: {alignments[0].on_same_strand()}'
                )
                if report:
                    # Extract zero-based coordinates of TIRs
                    ref_start = alignments[0].ref_coords().start
                    ref_end = alignments[0].ref_coords().end
                    # qry_coords are in correct left to right order, even if hit is on opposite strand.
                    qry_start = alignments[0].qry_coords().start
                    qry_end = alignments[0].qry_coords().end

                    # hit_length_ref appears to be the length of the entire aligned ref seq including gaps.
                    # hit_length_qry appears to be only the matching bases between the two TIRs.
                    # if alignments[0].hit_length_ref != alignments[0].hit_length_qry:
                    #    logging.warning(
                    #        f'Length of TIRs do not match: {alignments[0].hit_length_ref}bp vs {alignments[0].hit_length_qry}bp \nref_coords: {alignments[0].ref_coords()} \nqry_coords: {alignments[0].qry_coords()} \nYou may want to inspect the alignment:'
                    #    )
                    # else:
                    #    logging.info('Alignment of TIRs:')

                    hit_length_ref = (ref_end + 1) - ref_start
                    hit_length_qry = (qry_end + 1) - qry_start

                    if hit_length_ref != hit_length_qry:
                        logging.warning(
                            f'Length of TIRs do not match: {hit_length_ref}bp vs {hit_length_qry}bp \nref_coords: {alignments[0].ref_coords()} \nqry_coords: {alignments[0].qry_coords()} \nYou may want to inspect the alignment:'
                        )
                    else:
                        logging.info('Alignment of TIRs:')

                    # Extract the ref and qry sequences
                    ref_seq = rec[ref_start : ref_end + 1].seq
                    qry_seq = rec[qry_start : qry_end + 1].seq

                    # Reverse complement the qry sequence
                    qry_seq = qry_seq.reverse_complement()

                    # Perform global alignment using PairwiseAligner
                    aligner = PairwiseAligner(scoring='blastn')
                    pairwise_alignments = aligner.align(ref_seq, qry_seq)

                    # Print the first gapped alignment
                    if pairwise_alignments:
                        print(pairwise_alignments[0])

                    if report in ['split', 'external', 'all']:
                        if both:
                            # Yield both left and right TIRs
                            # Left TIR
                            leftSeg = rec[ref_start : ref_end + 1]
                            leftSeg.id = f'{leftSeg.id}_L_TIR'
                            leftSeg.name = leftSeg.id
                            leftSeg.description = f'[{rec.id} left TIR segment]'
                            logging.info(
                                f'Yielding left TIR segment: {leftSeg.id}, len: {len(leftSeg)}bp'
                            )
                            yield leftSeg

                            # Right TIR - reverse complement so it aligns with left TIR
                            rightSeg = rec[qry_start : qry_end + 1]
                            # Reverse complement the right TIR
                            rightSeg = rightSeg.reverse_complement(
                                id=f'{rec.id}_R_TIR',
                                name=f'{rec.id}_R_TIR',
                                description=f'[{rec.id} right TIR segment, reverse complemented]',
                            )
                            logging.info(
                                f'Yielding right TIR segment (reverse complemented): {rightSeg.id}, len: {len(rightSeg)}bp'
                            )
                            yield rightSeg
                        else:
                            # yield TIR slice - append "_TIR"
                            extSeg = rec[
                                ref_start : ref_end + 1
                            ]  # +1 to include end base
                            extSeg.id = f'{extSeg.id}_TIR'
                            extSeg.name = extSeg.id
                            extSeg.description = f'[{rec.id} TIR segment]'
                            logging.info(
                                f'Yielding TIR segment: {extSeg.id}, len: {len(extSeg)}bp'
                            )
                            yield extSeg

                    if report in ['split', 'internal', 'all']:
                        # yield internal slice - append "_I"
                        intSeg = rec[
                            ref_end + 1 : qry_start
                        ]  # no +1 as we want to exclude the base at qry_start
                        intSeg.id = f'{intSeg.id}_I'
                        intSeg.name = intSeg.id
                        intSeg.description = f'[{rec.id} internal segment]'
                        logging.info(
                            f'Yielding internal segment: {intSeg.id}, len: {len(intSeg)}bp'
                        )
                        yield intSeg
                    if report == 'all':
                        logging.info(
                            f'Yielding original element: {rec.id}, len: {len(rec)}bp'
                        )
                        yield rec
                if mites:
                    # Assemble TIRs into hypothetical MITEs
                    spacer = 0  # Spacer length between TIRs. Future option.
                    synMITE = (
                        rec[ref_start : ref_end + 1 + spacer]
                        + rec[qry_start : qry_end + 1]
                    )
                    synMITE.id = f'{synMITE.id}_synMITE'
                    synMITE.name = synMITE.id
                    synMITE.description = (
                        f'[Synthetic MITE constructed from {rec.id} TIRs]'
                    )
                    logging.info(f'Yielding synthetic MITE: {synMITE.id}')
                    yield synMITE
            else:
                # If alignment list empty after filtering, print alert and continue
                logging.info(f'No TIRs found for candidate element: {rec.id}')
    finally:
        # Log finished processing len(seen_ids) elements
        logging.info(f'Finished processing {len(seen_ids)} elements.')

        # Log number of elements with TIRs
        logging.info(f'Found TIRs in {found_TIRs_count} elements.')

        # Clean up the temporary directory if keeptemp is False
        if not keeptemp:
            shutil.rmtree(tempDir)
            logging.info(f'Temporary directory deleted: {tempDir}')
        else:
            logging.info(f'Temporary directory retained: {tempDir}')


def filterCoordsFileTIR(
    coordsFile: str, record: SeqRecord, minterm: int = 10, flankdist: int = 10
) -> List:
    """
    Filter alignment coordinates file for TIR candidates.

    Processes alignment results to identify Terminal Inverted Repeats (TIRs)
    based on position, orientation, and length criteria.

    Parameters
    ----------
    coordsFile : str
        Path to the coords file with alignment information.
    record : SeqRecord
        The sequence record being analyzed.
    minterm : int, optional
        Minimum length for a terminal repeat to be considered, by default 10.
    flankdist : int, optional
        Maximum distance from element boundaries for TIR candidates, by default 10.

    Returns
    -------
    list
        List of alignment objects that meet TIR criteria, sorted by
        internal segment size (largest first).

    Notes
    -----
    Filtering steps include:
    1. Removing self-alignments
    2. Keeping alignments meeting minimum length requirement
    3. Keeping alignments on opposite strands (characteristic of TIRs)
    4. Filtering for alignments near element boundaries
    5. Filtering out overlapping alignments
    6. Sorting by internal segment size
    """
    # Import coords file to iterator object
    file_reader = coords_file.reader(coordsFile)

    # Exclude hits to self. Also converts iterator output to stable list
    alignments = [hit for hit in file_reader if not hit.is_self_hit()]

    logging.debug(f'NON SELF ALIGNMENTS: {len(alignments)}')

    # Filter hits less than min length (Done internally for nucmer, not blastn.)
    alignments = [hit for hit in alignments if hit.hit_length_ref >= minterm]

    logging.debug(f'ALIGNMENTS >= minlen {minterm}bp: {len(alignments)}')

    # Filter for hits on same strand i.e. tandem repeats / LTRs
    alignments = [hit for hit in alignments if not hit.on_same_strand()]

    logging.debug(f'ALIGNMENTS ON OPPOSITE STRANDS: {len(alignments)}')

    # Filter for 5' repeats which begin within x bases of element start
    # hit.ref_start is zero-based, so we don't need to add 1
    # flankdist zero enforces that the hit starts at the beginning of the element
    alignments = [hit for hit in alignments if hit.ref_start <= flankdist]

    logging.debug(
        f'ALIGNMENTS within {flankdist}bp of element start: {len(alignments)}'
    )

    # Filter for 5' repeats with complementary hits within x bases of element end
    # Note: when hit is on opposite strand, hit.qry_start is the 3' end
    alignments = [
        hit for hit in alignments if (len(record) - (hit.qry_start + 1)) <= flankdist
    ]

    logging.debug(
        f'ALIGNMENTS with hit within {flankdist}bp of element end at {len(record)}bp: {len(alignments)}'
    )

    # Scrub overlapping ref / query segments
    alignments = [hit for hit in alignments if hit.ref_end < hit.qry_end]

    logging.debug(f'NON-OVERLAPPING ALIGNMENTS: {len(alignments)}')

    # Sort largest to smallest dist between end of ref (subject) and start
    # of query (hit)
    # x.qry_end - x.ref_end =
    # 5'end of right TIR - 3' end of left TIR = length of internal segment
    # TIR pair with largest internal segment (outermost TIRs) is first in list.
    alignments = sorted(alignments, key=lambda x: (x.qry_end - x.ref_end), reverse=True)

    return alignments


def getLTRs(
    fasta_file: str,
    flankdist: int = 10,
    minid: float = 80,
    minterm: int = 10,
    minseed: int = 5,
    diagfactor: float = 0.3,
    report: str = 'split',
    temp: Optional[str] = None,
    keeptemp: bool = False,
    alignTool: str = 'nucmer',
    verbose: bool = True,
    both: bool = False,
) -> Generator[SeqRecord, None, None]:
    """
    Align elements to self and attempt to identify LTRs.

    Processes sequences to identify Long Terminal Repeats (LTRs)
    by performing self-alignment and filtering results.

    Parameters
    ----------
    fasta_file : str
        Path to the multifasta file containing sequence records.
    flankdist : int, optional
        Maximum distance from element start for LTR candidates, by default 10.
    minid : float, optional
        Minimum identity between terminal repeat pairs, by default 80.
    minterm : int, optional
        Minimum length for a terminal repeat to be considered, by default 10.
    minseed : int, optional
        Minimum seed length for nucmer, by default 5.
    diagfactor : float, optional
        Diagonal factor for nucmer, by default 0.3.
    report : str, optional
        Reporting mode for LTRs ('split', 'external', 'internal', 'all'), by default 'split'.
    temp : str, optional
        Path to the temporary directory, by default None.
    keeptemp : bool, optional
        Whether to keep the temporary directory after processing, by default False.
    alignTool : str, optional
        Alignment tool to use ('nucmer' or 'blastn'), by default 'nucmer'.
    verbose : bool, optional
        Whether to print verbose output, by default True.
    both : bool, optional
        Whether to report both left and right terminal repeats, by default False.

    Yields
    ------
    SeqRecord
        Segments of the sequence based on the reporting mode:
        - 'split': LTRs and internal regions separately
        - 'external': Only LTRs
        - 'internal': Only internal regions
        - 'all': Original sequences plus all segments

    Notes
    -----
    When both=True and report is in ['split', 'external', 'all'], both left and right
    LTRs will be yielded with suffixes '_L_LTR' and '_R_LTR'. Unlike TIRs, the right
    LTR will NOT be reverse complemented as it is already in the same orientation as
    the left LTR.
    """
    # Set temp directory to cwd if none is provided
    if not temp:
        temp = os.getcwd()

    # Create a unique temporary directory
    tempDir = tempfile.mkdtemp(prefix=f'tsplit_LTR_temp_{alignTool}_', dir=temp)
    logging.debug(f'Temporary directory created: {tempDir}')

    # Set to store seen IDs
    seen_ids = set()

    # Count of elements with LTRs
    found_LTRs_count = 0

    try:
        # Iterate over each record in the fasta file
        for rec in SeqIO.parse(fasta_file, 'fasta'):
            # Log the record name and length
            logging.info(f'Processing record: {rec.id}, Length: {len(rec)}bp')

            # Check for duplicate IDs
            if rec.id in seen_ids:
                logging.error(f'Duplicate record ID found: {rec.id}')
                raise ValueError(f'Duplicate record ID found: {rec.id}')

            seen_ids.add(rec.id)

            # Create temp paths for single element fasta and alignment coords
            tempFasta = os.path.join(tempDir, cleanID(rec.id) + '.fasta')
            tempCoords = os.path.join(tempDir, f'{cleanID(rec.id)}_{alignTool}.coords')

            # Write current element to single fasta
            with open(tempFasta, 'w') as f:
                SeqIO.write(rec, f, 'fasta')

            # Align to self with nucmer
            if alignTool == 'nucmer':
                # Compose Nucmer script for current element vs self
                runner = nucmer.Runner(
                    tempFasta,
                    tempFasta,
                    tempCoords,
                    min_id=minid,
                    min_length=minseed,
                    diagfactor=diagfactor,
                    mincluster=minterm,
                    breaklen=200,
                    maxmatch=True,
                    simplify=False,
                )
                # Execute nucmer
                runner.run()
            elif alignTool == 'blastn':
                # Alternatively, use blastn as search tool and write nucmer.coords-like output.
                cmd = makeBlast(seq=tempFasta, outfile=tempCoords, pid=minid)
                run_cmd(cmd, verbose=verbose, workingDir=tempDir)

            alignments = filterCoordsFileLTR(
                tempCoords, rec, minterm=minterm, flankdist=flankdist
            )

            # If alignments exist after filtering report features using alignment pair with largest
            # internal segment i.e. first element in sorted list.
            if alignments:
                # Increment the count of elements with LTRs
                found_LTRs_count += 1

                logging.info(
                    f'Found {len(alignments)} alignments for candidate element: {rec.id}'
                )

                [print(x) for x in alignments]

                logging.info('Selecting candidate LTRs with largest internal segment:')

                print(alignments[0])

                logging.debug(
                    f'Identified coords in top candidate: \nref_start: {alignments[0].ref_start} \nref_end: {alignments[0].ref_end} \nqry_start: {alignments[0].qry_start} \nqry_end: {alignments[0].qry_end} \nref_coords: {alignments[0].ref_coords()} \nqry_coords: {alignments[0].qry_coords()} \nref_length: {alignments[0].ref_length} \nqry_length: {alignments[0].qry_length} \nframe: {alignments[0].frame} \npercent_identity: {alignments[0].percent_identity} \nhit_length_ref: {alignments[0].hit_length_ref} \nhit_length_qry: {alignments[0].hit_length_qry} \non_same_strand: {alignments[0].on_same_strand()}'
                )

                if report:
                    # Extract zero-based coordinates of LTRs
                    ref_start = alignments[0].ref_coords().start
                    ref_end = alignments[0].ref_coords().end
                    qry_start = alignments[0].qry_coords().start
                    qry_end = alignments[0].qry_coords().end

                    if alignments[0].hit_length_ref != alignments[0].hit_length_qry:
                        logging.warning(
                            f'Length of LTRs do not match: {alignments[0].hit_length_ref}bp vs {alignments[0].hit_length_qry}bp \nref_coords: {alignments[0].ref_coords()} \nqry_coords: {alignments[0].qry_coords()} \nYou may want to inspect the alignment:'
                        )
                    else:
                        logging.info('Alignment of LTRs:')

                    # Extract the ref and qry sequences
                    ref_seq = rec[ref_start : ref_end + 1].seq
                    qry_seq = rec[qry_start : qry_end + 1].seq

                    # Perform global alignment using PairwiseAligner
                    aligner = PairwiseAligner(scoring='blastn')
                    pairwise_alignments = aligner.align(ref_seq, qry_seq)

                    # Print the first gapped alignment
                    if pairwise_alignments:
                        print(pairwise_alignments[0])

                    if report in ['split', 'external', 'all']:
                        if both:
                            # Yield both left and right LTRs
                            # Left LTR
                            leftSeg = rec[ref_start : ref_end + 1]
                            leftSeg.id = f'{leftSeg.id}_L_LTR'
                            leftSeg.name = leftSeg.id
                            leftSeg.description = f'[{rec.id} left LTR segment]'
                            logging.info(
                                f'Yielding left LTR segment: {leftSeg.id}, len: {len(leftSeg)}bp'
                            )
                            yield leftSeg

                            # Right LTR - keep in same orientation (do NOT reverse complement)
                            rightSeg = rec[qry_start : qry_end + 1]
                            rightSeg.id = f'{rightSeg.id}_R_LTR'
                            rightSeg.name = rightSeg.id
                            rightSeg.description = f'[{rec.id} right LTR segment]'
                            logging.info(
                                f'Yielding right LTR segment: {rightSeg.id}, len: {len(rightSeg)}bp'
                            )
                            yield rightSeg
                        else:
                            # yield LTR slice - append "_LTR"
                            extSeg = rec[
                                ref_start : ref_end + 1
                            ]  # +1 to include end base
                            extSeg.id = f'{extSeg.id}_LTR'
                            extSeg.name = extSeg.id
                            extSeg.description = f'[{rec.id} LTR segment]'
                            logging.info(
                                f'Yielding LTR segment: {extSeg.id}, len: {len(extSeg)}bp'
                            )
                            yield extSeg

                    if report in ['split', 'internal', 'all']:
                        # yield internal slice - append "_I"
                        intSeg = rec[
                            ref_end + 1 : qry_start
                        ]  # no +1 as we want to exclude the base at qry_start
                        intSeg.id = f'{intSeg.id}_I'
                        intSeg.name = intSeg.id
                        intSeg.description = f'[{rec.id} internal segment]'
                        logging.info(
                            f'Yielding internal segment: {intSeg.id}, len: {len(intSeg)}bp'
                        )
                        yield intSeg

                    if report == 'all':
                        # yield original element
                        logging.info(
                            f'Yielding original element: {rec.id}, len: {len(rec)}bp'
                        )
                        yield rec
            else:
                # If alignment list is empty after filtering, print alert and continue.
                logging.info(f'No LTRs found for candidate element: {rec.id}')
    finally:
        # Log finished processing len(seen_ids) elements
        logging.info(f'Finished processing {len(seen_ids)} elements.')

        # Log number of elements with LTRs
        logging.info(f'Found LTRs in {found_LTRs_count} elements.')

        # Clean up the temporary directory if keeptemp is False
        if not keeptemp:
            shutil.rmtree(tempDir)
            logging.info(f'Cleaning up temporary directory: {tempDir}')
        else:
            logging.info(f'Temporary directory retained: {tempDir}')


def filterCoordsFileLTR(
    coordsFile: str, record: SeqRecord, minterm: int = 10, flankdist: int = 10
) -> List:
    """
    Filter alignment coordinates file for LTR candidates.

    Processes alignment results to identify Long Terminal Repeats (LTRs)
    based on position, orientation, and length criteria.

    Parameters
    ----------
    coordsFile : str
        Path to the coords file with alignment information.
    record : SeqRecord
        The sequence record being analyzed.
    minterm : int, optional
        Minimum length for a terminal repeat to be considered, by default 10.
    flankdist : int, optional
        Maximum distance from element boundaries for LTR candidates, by default 10.

    Returns
    -------
    list
        List of alignment objects that meet LTR criteria, sorted by
        internal segment size (largest first).

    Notes
    -----
    Filtering steps include:
    1. Removing self-alignments
    2. Keeping alignments meeting minimum length requirement
    3. Keeping alignments on the same strand (characteristic of LTRs)
    4. Filtering for alignments near element boundaries
    5. Filtering out overlapping alignments
    6. Sorting by internal segment size
    """
    # Import coords file to iterator object
    file_reader = coords_file.reader(coordsFile)

    # Exclude hits to self. Also converts iterator output to stable list
    alignments = [hit for hit in file_reader if not hit.is_self_hit()]

    logging.debug(f'NON SELF ALIGNMENTS: {len(alignments)}')

    # Filter hits less than min length (Done internally for nucmer, not blastn.)
    alignments = [hit for hit in alignments if hit.hit_length_ref >= minterm]

    logging.debug(f'ALIGNMENTS >= minlen {minterm}bp: {len(alignments)}')

    # Filter for hits on same strand i.e. tandem repeats / LTRs
    alignments = [hit for hit in alignments if hit.on_same_strand()]

    logging.debug(f'ALIGNMENTS ON SAME STRAND: {len(alignments)}')

    # Filter for 5' repeats which begin within x bases of element start
    # hit.ref_start is zero-based, so we don't need to add 1
    # flankdist zero enforces that the hit starts at the beginning of the element
    alignments = [hit for hit in alignments if hit.ref_start <= flankdist]

    logging.debug(
        f'ALIGNMENTS within {flankdist}bp of element start: {len(alignments)}'
    )

    # Filter for 5' repeats whose 3' match ends within x bases of element end
    alignments = [
        hit for hit in alignments if len(record) - hit.qry_end + 1 <= flankdist
    ]

    logging.debug(
        f'ALIGNMENTS with hit within {flankdist}bp of element end at {len(record)}bp: {len(alignments)}'
    )

    # Keep non-overlappying ref / query segments
    alignments = [hit for hit in alignments if hit.ref_end < hit.qry_start]

    logging.debug(f'NON-OVERLAPPING ALIGNMENTS: {len(alignments)}')

    # Sort largest to smallest dist between end of ref (subject) and start of query (hit)
    # x.qry_start (3') - x.ref_end (5') = Length of internal segment
    alignments = sorted(
        alignments, key=lambda x: (x.qry_start - x.ref_end), reverse=True
    )

    return alignments


def getTIRs_with_data(
    fasta_file: str,
    flankdist: int = 10,
    minid: float = 80,
    minterm: int = 10,
    minseed: int = 5,
    diagfactor: float = 0.3,
    mites: bool = False,
    report: str = 'split',
    temp: Optional[str] = None,
    keeptemp: bool = False,
    alignTool: str = 'nucmer',
    verbose: bool = True,
    both: bool = False,
    collect_alignments: bool = False,
    collect_features: bool = False,
) -> tuple:
    """
    Collect alignment and feature data while extracting TIRs.

    This function calls getTIRs and captures alignment data for PAF output
    and feature data for GFF3 output while preserving the original behavior.

    Parameters
    ----------
    All parameters are the same as getTIRs, plus:
    collect_alignments : bool, optional
        Whether to collect alignment data for PAF output, by default False.
    collect_features : bool, optional
        Whether to collect feature data for GFF3 output, by default False.

    Returns
    -------
    tuple
        (segments_generator, alignment_data, feature_data)
        - segments_generator: Generator yielding SeqRecords (same as getTIRs)
        - alignment_data: List of dicts with alignment info (empty if collect_alignments=False)
        - feature_data: List of dicts with feature info (empty if collect_features=False)
    """
    # Storage for alignment and feature data
    alignment_data = []
    feature_data = []

    # If not collecting data, just return the original generator
    if not collect_alignments and not collect_features:
        segments = getTIRs(
            fasta_file,
            flankdist=flankdist,
            minid=minid,
            minterm=minterm,
            minseed=minseed,
            diagfactor=diagfactor,
            mites=mites,
            report=report,
            temp=temp,
            keeptemp=keeptemp,
            alignTool=alignTool,
            verbose=verbose,
            both=both,
        )
        return segments, alignment_data, feature_data

    # We need to process the file ourselves to collect data
    # Set temp directory to cwd if none is provided
    if not temp:
        temp = os.getcwd()

    # Create a unique temporary directory
    tempDir = tempfile.mkdtemp(prefix=f'tsplit_TIR_temp_{alignTool}_', dir=temp)
    logging.debug(f'Temporary directory created: {tempDir}')

    seen_ids = set()
    segments_list = []

    try:
        # Iterate over each record in the fasta file
        for rec in SeqIO.parse(fasta_file, 'fasta'):
            logging.info(
                f'Processing record {len(seen_ids) + 1}: Name: {rec.id}, Length: {len(rec)}bp'
            )

            # Check for duplicate IDs
            if rec.id in seen_ids:
                logging.error(f'Duplicate record ID found: {rec.id}')
                raise ValueError(f'Duplicate record ID found: {rec.id}')

            seen_ids.add(rec.id)

            # Create temp paths for single element fasta and alignment coords
            tempFasta = os.path.join(tempDir, cleanID(rec.id) + '.fasta')
            tempCoords = os.path.join(
                tempDir, cleanID(rec.id) + '_' + alignTool + '.coords'
            )

            # Write current element to single fasta
            with open(tempFasta, 'w') as f:
                SeqIO.write(rec, f, 'fasta')

            # Align to self
            if alignTool == 'nucmer':
                runner = nucmer.Runner(
                    tempFasta,
                    tempFasta,
                    tempCoords,
                    min_id=minid,
                    min_length=minseed,
                    diagfactor=diagfactor,
                    mincluster=minterm,
                    breaklen=200,
                    maxmatch=True,
                    simplify=False,
                )
                runner.run()
            elif alignTool == 'blastn':
                cmd = makeBlast(seq=tempFasta, outfile=tempCoords, pid=minid)
                run_cmd(cmd, verbose=verbose, workingDir=tempDir)

            alignments = filterCoordsFileTIR(
                tempCoords, rec, minterm=minterm, flankdist=flankdist
            )

            # Collect ALL alignment data for PAF if requested (before filtering)
            if collect_alignments:
                # Read all alignments from coords file (including self-hits for complete PAF)
                file_reader = coords_file.reader(tempCoords)
                all_alignments = list(file_reader)

                # Add all alignments to PAF output
                for aln in all_alignments:
                    # Get coordinates from alignment
                    aln_ref_start = aln.ref_coords().start
                    aln_ref_end = aln.ref_coords().end
                    aln_qry_start = aln.qry_coords().start
                    aln_qry_end = aln.qry_coords().end

                    # Calculate strand
                    strand = '-' if not aln.on_same_strand() else '+'

                    # Calculate number of matches from percent identity
                    num_matches = int(
                        aln.hit_length_ref * aln.percent_identity / 100.0
                    )

                    alignment_data.append(
                        {
                            'qry_name': rec.id,
                            'qry_length': len(rec),
                            'qry_start': aln_qry_start,
                            'qry_end': aln_qry_end + 1,  # PAF uses end-exclusive
                            'strand': strand,
                            'ref_name': rec.id,
                            'ref_length': len(rec),
                            'ref_start': aln_ref_start,
                            'ref_end': aln_ref_end + 1,  # PAF uses end-exclusive
                            'num_matches': num_matches,
                            'aln_block_length': aln.hit_length_ref,
                            'mapping_quality': 255,  # Not available, use 255 (missing)
                        }
                    )

            # If alignments exist, collect data and generate segments
            if alignments:
                # Get coordinates from best alignment (first in sorted list)
                best_aln = alignments[0]
                ref_start = best_aln.ref_coords().start
                ref_end = best_aln.ref_coords().end
                qry_start = best_aln.qry_coords().start
                qry_end = best_aln.qry_coords().end

                # Collect feature data for GFF3 if requested
                if collect_features:
                    # Create GFF3 features for TIRs
                    # Convert to 1-based coordinates for GFF3
                    left_tir_start = ref_start + 1
                    left_tir_end = ref_end + 1
                    right_tir_start = qry_start + 1
                    right_tir_end = qry_end + 1

                    # Left TIR feature
                    feature_data.append(
                        {
                            'seqid': rec.id,
                            'type': 'TIR',
                            'start': left_tir_start,
                            'end': left_tir_end,
                            'score': f'{best_aln.percent_identity:.2f}',
                            'strand': '+',
                            'phase': '.',
                            'attributes': f'ID={rec.id}_TIR_L;Name=Left_TIR;length={ref_end - ref_start + 1};identity={best_aln.percent_identity:.2f}',
                        }
                    )

                    # Right TIR feature (on opposite strand for TIRs)
                    feature_data.append(
                        {
                            'seqid': rec.id,
                            'type': 'TIR',
                            'start': right_tir_start,
                            'end': right_tir_end,
                            'score': f'{best_aln.percent_identity:.2f}',
                            'strand': '-',
                            'phase': '.',
                            'attributes': f'ID={rec.id}_TIR_R;Name=Right_TIR;length={qry_end - qry_start + 1};identity={best_aln.percent_identity:.2f}',
                        }
                    )

                # Generate segments using original logic
                if report:
                    if report in ['split', 'external', 'all']:
                        if both:
                            leftSeg = rec[ref_start : ref_end + 1]
                            leftSeg.id = f'{leftSeg.id}_L_TIR'
                            leftSeg.name = leftSeg.id
                            leftSeg.description = f'[{rec.id} left TIR segment]'
                            segments_list.append(leftSeg)

                            rightSeg = rec[qry_start : qry_end + 1]
                            rightSeg = rightSeg.reverse_complement(
                                id=f'{rec.id}_R_TIR',
                                name=f'{rec.id}_R_TIR',
                                description=f'[{rec.id} right TIR segment, reverse complemented]',
                            )
                            segments_list.append(rightSeg)
                        else:
                            extSeg = rec[ref_start : ref_end + 1]
                            extSeg.id = f'{extSeg.id}_TIR'
                            extSeg.name = extSeg.id
                            extSeg.description = f'[{rec.id} TIR segment]'
                            segments_list.append(extSeg)

                    if report in ['split', 'internal', 'all']:
                        intSeg = rec[ref_end + 1 : qry_start]
                        intSeg.id = f'{intSeg.id}_I'
                        intSeg.name = intSeg.id
                        intSeg.description = f'[{rec.id} internal segment]'
                        segments_list.append(intSeg)

                    if report == 'all':
                        segments_list.append(rec)

                if mites:
                    spacer = 0
                    synMITE = (
                        rec[ref_start : ref_end + 1 + spacer]
                        + rec[qry_start : qry_end + 1]
                    )
                    synMITE.id = f'{synMITE.id}_synMITE'
                    synMITE.name = synMITE.id
                    synMITE.description = (
                        f'[Synthetic MITE constructed from {rec.id} TIRs]'
                    )
                    segments_list.append(synMITE)

    finally:
        # Clean up the temporary directory if keeptemp is False
        if not keeptemp:
            shutil.rmtree(tempDir)
            logging.info(f'Temporary directory deleted: {tempDir}')
        else:
            logging.info(f'Temporary directory retained: {tempDir}')

    # Return segments as a generator
    return iter(segments_list), alignment_data, feature_data


def getLTRs_with_data(
    fasta_file: str,
    flankdist: int = 10,
    minid: float = 80,
    minterm: int = 10,
    minseed: int = 5,
    diagfactor: float = 0.3,
    report: str = 'split',
    temp: Optional[str] = None,
    keeptemp: bool = False,
    alignTool: str = 'nucmer',
    verbose: bool = True,
    both: bool = False,
    collect_alignments: bool = False,
    collect_features: bool = False,
) -> tuple:
    """
    Collect alignment and feature data while extracting LTRs.

    This function calls getLTRs and captures alignment data for PAF output
    and feature data for GFF3 output while preserving the original behavior.

    Parameters
    ----------
    All parameters are the same as getLTRs, plus:
    collect_alignments : bool, optional
        Whether to collect alignment data for PAF output, by default False.
    collect_features : bool, optional
        Whether to collect feature data for GFF3 output, by default False.

    Returns
    -------
    tuple
        (segments_generator, alignment_data, feature_data)
        - segments_generator: Generator yielding SeqRecords (same as getLTRs)
        - alignment_data: List of dicts with alignment info (empty if collect_alignments=False)
        - feature_data: List of dicts with feature info (empty if collect_features=False)
    """
    # Storage for alignment and feature data
    alignment_data = []
    feature_data = []

    # If not collecting data, just return the original generator
    if not collect_alignments and not collect_features:
        segments = getLTRs(
            fasta_file,
            flankdist=flankdist,
            minid=minid,
            minterm=minterm,
            minseed=minseed,
            diagfactor=diagfactor,
            report=report,
            temp=temp,
            keeptemp=keeptemp,
            alignTool=alignTool,
            verbose=verbose,
            both=both,
        )
        return segments, alignment_data, feature_data

    # We need to process the file ourselves to collect data
    # Set temp directory to cwd if none is provided
    if not temp:
        temp = os.getcwd()

    # Create a unique temporary directory
    tempDir = tempfile.mkdtemp(prefix=f'tsplit_LTR_temp_{alignTool}_', dir=temp)
    logging.debug(f'Temporary directory created: {tempDir}')

    seen_ids = set()
    segments_list = []

    try:
        # Iterate over each record in the fasta file
        for rec in SeqIO.parse(fasta_file, 'fasta'):
            logging.info(
                f'Processing record {len(seen_ids) + 1}: Name: {rec.id}, Length: {len(rec)}bp'
            )

            # Check for duplicate IDs
            if rec.id in seen_ids:
                logging.error(f'Duplicate record ID found: {rec.id}')
                raise ValueError(f'Duplicate record ID found: {rec.id}')

            seen_ids.add(rec.id)

            # Create temp paths for single element fasta and alignment coords
            tempFasta = os.path.join(tempDir, cleanID(rec.id) + '.fasta')
            tempCoords = os.path.join(
                tempDir, cleanID(rec.id) + '_' + alignTool + '.coords'
            )

            # Write current element to single fasta
            with open(tempFasta, 'w') as f:
                SeqIO.write(rec, f, 'fasta')

            # Align to self
            if alignTool == 'nucmer':
                runner = nucmer.Runner(
                    tempFasta,
                    tempFasta,
                    tempCoords,
                    min_id=minid,
                    min_length=minseed,
                    diagfactor=diagfactor,
                    mincluster=minterm,
                    breaklen=200,
                    maxmatch=True,
                    simplify=False,
                )
                runner.run()
            elif alignTool == 'blastn':
                cmd = makeBlast(seq=tempFasta, outfile=tempCoords, pid=minid)
                run_cmd(cmd, verbose=verbose, workingDir=tempDir)

            alignments = filterCoordsFileLTR(
                tempCoords, rec, minterm=minterm, flankdist=flankdist
            )

            # Collect ALL alignment data for PAF if requested (before filtering)
            if collect_alignments:
                # Read all alignments from coords file (including self-hits for complete PAF)
                file_reader = coords_file.reader(tempCoords)
                all_alignments = list(file_reader)

                # Add all alignments to PAF output
                for aln in all_alignments:
                    # Get coordinates from alignment
                    aln_ref_start = aln.ref_coords().start
                    aln_ref_end = aln.ref_coords().end
                    aln_qry_start = aln.qry_coords().start
                    aln_qry_end = aln.qry_coords().end

                    # Calculate strand
                    strand = '+' if aln.on_same_strand() else '-'

                    # Calculate number of matches from percent identity
                    num_matches = int(
                        aln.hit_length_ref * aln.percent_identity / 100.0
                    )

                    alignment_data.append(
                        {
                            'qry_name': rec.id,
                            'qry_length': len(rec),
                            'qry_start': aln_qry_start,
                            'qry_end': aln_qry_end + 1,  # PAF uses end-exclusive
                            'strand': strand,
                            'ref_name': rec.id,
                            'ref_length': len(rec),
                            'ref_start': aln_ref_start,
                            'ref_end': aln_ref_end + 1,  # PAF uses end-exclusive
                            'num_matches': num_matches,
                            'aln_block_length': aln.hit_length_ref,
                            'mapping_quality': 255,  # Not available, use 255 (missing)
                        }
                    )

            # If alignments exist, collect data and generate segments
            if alignments:
                # Get coordinates from best alignment (first in sorted list)
                best_aln = alignments[0]
                ref_start = best_aln.ref_coords().start
                ref_end = best_aln.ref_coords().end
                qry_start = best_aln.qry_coords().start
                qry_end = best_aln.qry_coords().end

                # Collect feature data for GFF3 if requested
                if collect_features:
                    # Create GFF3 features for LTRs
                    # Convert to 1-based coordinates for GFF3
                    left_ltr_start = ref_start + 1
                    left_ltr_end = ref_end + 1
                    right_ltr_start = qry_start + 1
                    right_ltr_end = qry_end + 1

                    # Left LTR feature
                    feature_data.append(
                        {
                            'seqid': rec.id,
                            'type': 'LTR',
                            'start': left_ltr_start,
                            'end': left_ltr_end,
                            'score': f'{best_aln.percent_identity:.2f}',
                            'strand': '+',
                            'phase': '.',
                            'attributes': f'ID={rec.id}_LTR_L;Name=Left_LTR;length={ref_end - ref_start + 1};identity={best_aln.percent_identity:.2f}',
                        }
                    )

                    # Right LTR feature (on same strand for LTRs)
                    feature_data.append(
                        {
                            'seqid': rec.id,
                            'type': 'LTR',
                            'start': right_ltr_start,
                            'end': right_ltr_end,
                            'score': f'{best_aln.percent_identity:.2f}',
                            'strand': '+',
                            'phase': '.',
                            'attributes': f'ID={rec.id}_LTR_R;Name=Right_LTR;length={qry_end - qry_start + 1};identity={best_aln.percent_identity:.2f}',
                        }
                    )

                # Generate segments using original logic
                if report:
                    if report in ['split', 'external', 'all']:
                        if both:
                            leftSeg = rec[ref_start : ref_end + 1]
                            leftSeg.id = f'{leftSeg.id}_L_LTR'
                            leftSeg.name = leftSeg.id
                            leftSeg.description = f'[{rec.id} left LTR segment]'
                            segments_list.append(leftSeg)

                            rightSeg = rec[qry_start : qry_end + 1]
                            rightSeg.id = f'{rec.id}_R_LTR'
                            rightSeg.name = rightSeg.id
                            rightSeg.description = f'[{rec.id} right LTR segment]'
                            segments_list.append(rightSeg)
                        else:
                            extSeg = rec[ref_start : ref_end + 1]
                            extSeg.id = f'{extSeg.id}_LTR'
                            extSeg.name = extSeg.id
                            extSeg.description = f'[{rec.id} LTR segment]'
                            segments_list.append(extSeg)

                    if report in ['split', 'internal', 'all']:
                        intSeg = rec[ref_end + 1 : qry_start]
                        intSeg.id = f'{intSeg.id}_I'
                        intSeg.name = intSeg.id
                        intSeg.description = f'[{rec.id} internal segment]'
                        segments_list.append(intSeg)

                    if report == 'all':
                        segments_list.append(rec)

    finally:
        # Clean up the temporary directory if keeptemp is False
        if not keeptemp:
            shutil.rmtree(tempDir)
            logging.info(f'Temporary directory deleted: {tempDir}')
        else:
            logging.info(f'Temporary directory retained: {tempDir}')

    # Return segments as a generator
    return iter(segments_list), alignment_data, feature_data

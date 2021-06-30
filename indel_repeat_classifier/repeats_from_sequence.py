import csv
import gzip
import logging
from pyfaidx import Fasta

LOGGER = logging.getLogger("RepeatFinder")
LOGGER.setLevel(logging.INFO)
formatter = logging.Formatter(
    '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel(LOGGER.level)
ch.setFormatter(formatter)
LOGGER.addHandler(ch)


def output_repeat(writer, seq, seqname, rpt_unit, start, end, flanks=10):
    if 'N' in rpt_unit:
        return 0
    end_pad = 0
    rpt_unit_length = len(rpt_unit)
    for i in range(rpt_unit_length - 1):
        tail = end + i
        if tail >= len(seq):
            break
        if seq[tail] == seq[tail - rpt_unit_length]:
            # e.g. if AG repeat with a final A: "nnnAGAGAnnn"
            end_pad += 1
        else:
            break
    end += end_pad
    f_start = max(0, start - flanks)
    flank_seq = seq[f_start:start].lower() + \
        seq[start:end].upper() + \
        seq[end:end + flanks].lower()
    rpt_length = end - start
    clps_len = 2 * rpt_unit_length
    if end_pad:  # we added partial repeat to end
        clps_len += end_pad
    clps_seq = seq[f_start:start].lower() + \
                seq[start:start + clps_len].upper() + \
                seq[end:end + flanks].lower()
    writer.writerow({
        'seqname': seqname,
        'start': start,
        'end': end,
        'rpt_type': 'Perfect',
        'rpt_length': rpt_length,
        'rpt_unit': rpt_unit,
        'collapsed_seq': clps_seq,
        'flank_seq': flank_seq
    })
    return 1


def update_progress(n, seqname, pos, progress_interval):
    if n % progress_interval == 0 and n != 0:
        LOGGER.info("Parsed {:,} repeats, at pos {}:{}".format(
            n, seqname, pos))


def find_repeats_in_string(seq,
                           seqname,
                           rpt_unit_length,
                           writer,
                           repeats_processed=0,
                           flanks=10,
                           progress_interval=1_000_000):
    window = rpt_unit_length * 2
    prev_rpt = None
    i = 0
    while i < len(seq) - rpt_unit_length:
        this_window = seq[i:i + window]
        if this_window[:rpt_unit_length] == this_window[rpt_unit_length:]:
            this_rpt = this_window[:rpt_unit_length]
            if prev_rpt == this_rpt:
                rpt_end = i + window
            else:
                if prev_rpt is not None:
                    repeats_processed += output_repeat(writer, seq, seqname,
                                                       prev_rpt, rpt_start,
                                                       rpt_end, flanks)
                    update_progress(repeats_processed, seqname, i,
                                    progress_interval)
                prev_rpt = this_rpt
                rpt_start = i
                rpt_end = i + window
            i += rpt_unit_length
        else:
            if prev_rpt is not None:
                repeats_processed += output_repeat(writer, seq, seqname,
                                                   prev_rpt, rpt_start,
                                                   rpt_end, flanks)
                update_progress(repeats_processed, seqname, i,
                                progress_interval)
            else:
                # check for microhomology
                mh_len = 0
                mh = None
                for j in range(1, rpt_unit_length):
                    x = this_window[:j]
                    y = this_window[rpt_unit_length:rpt_unit_length + j]
                    if x == y:
                        mh_len = j + rpt_unit_length
                        mh = this_window[:rpt_unit_length]
                    else:
                        break
                if mh is not None and 'N' not in mh:
                    mh_end = i + mh_len
                    f_start = max(0, i - flanks)
                    flank_seq = seq[f_start:i].lower() + \
                        seq[i:mh_end].upper() + \
                        seq[mh_end:mh_end + flanks].lower()
                    writer.writerow({
                        'seqname': seqname,
                        'start': i,
                        'end': i + mh_len,
                        'rpt_type': 'Imperfect',
                        'rpt_length': mh_len,
                        'rpt_unit': mh,
                        'collapsed_seq': flank_seq,
                        'flank_seq': flank_seq
                    })
                    repeats_processed += 1
            prev_rpt = None
            rpt_start = None
            rpt_end = None
            i += 1
    if prev_rpt is not None:
        repeats_processed += output_repeat(writer, seq, seqname, prev_rpt,
                                           rpt_start, rpt_end, flanks)
    LOGGER.info("Finished parsing {:,} repeats".format(repeats_processed))
    return repeats_processed


def find_repeats_in_fasta(fasta_path,
                          out_path,
                          rpt_unit_length,
                          flanks=10,
                          seqnames=[],
                          blacklist=[],
                          progress_interval=1_000_000):
    '''
    Output repeats from fasta file to CSV.

    Args:
            fasta_path:
                    Path to fasta file to process

            out_path:
                    Path to CSV output. Will be gzip compressed.

            rpt_unit_length:
                    Find repeats of this unit length. Set to 2 to find
                    dinucleotide repeats (e.g. AGAG...), set to 3 to find
                    trinucleotide repeats (e.g. AGTAGT...) etc.

            flanks:
                    Output this many bp of sequence either side of each repeat.

            sequence:
                    Process these sequences only.

            blacklist:
                    Skip these sequence names from fasta file.

            progress_interval:
                    Log progress once this many repeats have been processed.
    '''
    fasta = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
    if not out_path.endswith('.gz'):
        out_path += '.gz'
    with gzip.open(out_path, 'wt', newline='') as csvfile:
        n = 0
        fieldnames = [
            'seqname', 'start', 'end', 'rpt_type', 'rpt_length', 'rpt_unit',
            'collapsed_seq', 'flank_seq'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        if not seqnames:
            seqnames = fasta.keys()
        for seq in seqnames:
            if blacklist and seq in blacklist:
                continue
            LOGGER.info("Processing sequence {} {}bp repeats".format(
                seq, rpt_unit_length))
            cseq = str(fasta[seq])
            n += find_repeats_in_string(seq=cseq,
                                        seqname=seq,
                                        rpt_unit_length=rpt_unit_length,
                                        flanks=flanks,
                                        writer=writer,
                                        repeats_processed=n,
                                        progress_interval=progress_interval)

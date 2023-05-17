import argparse

from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
import numpy as np


def rna_to_dna(seq):
    return seq.upper().replace("U", "T")


def revcomp(seq):
    return "".join(reversed(seq.translate(str.maketrans("ATGC", "TACG"))))


def parse_range(range_txt):
    numbers = list()
    chunks = range_txt.replace(" ", "").split(",")
    for chunk in chunks:
        if "-" in chunk:
            start, end = chunk.split("-")
            for i in range(int(start), int(end) + 1):
                assert i not in numbers
                numbers.append(i)
        else:
            numbers.append(int(chunk))
    return sorted(numbers)

assert parse_range("1") == [1]
assert parse_range("1-6") == [1, 2, 3, 4, 5, 6]
assert parse_range("19-22, 10 - 12") == [10, 11, 12, 19, 20, 21, 22]


def read_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) != 1:
        raise ValueError("Currently only fasta files with one record are supported.")
    seq = str(records[0].seq)
    return seq


def grow_primer(template, start_pos, primer_type, growth_direction, min_tm, min_len=20):
    primer_tm = None
    primer_len = min_len
    primer_seq = None
    while primer_tm is None or primer_tm < min_tm:
        if growth_direction == "dn":
            # grow primer in downstream direction relative to template
            lower, upper = start_pos, primer_len + start_pos - 1
        elif growth_direction == "up":
            # grow primer in upstream direction relative to template
            lower, upper = start_pos - primer_len + 1, start_pos
        else:
            raise ValueError(growth_direction)
        if lower < 1 or upper > len(template):
            raise ValueError("Could not find suitable primer.")
        if primer_type == "fwd":
            primer_seq = template[lower - 1: upper]
        elif primer_type == "rev":
            primer_seq = revcomp(template[lower - 1: upper])
        else:
            raise ValueError(primer_type)
        primer_tm = mt.Tm_NN(primer_seq, nn_table=mt.DNA_NN2, saltcorr=7)
        primer_len += 1
    return primer_seq


def olaygo(input_file, output_file, tile_len_min, tile_len_max, primer_span_max, primer_tm_min=60.0, oligo_range=None, gap=2):
    seq = rna_to_dna(read_fasta(input_file))
    if oligo_range is None:
        target_pos = np.arange(1, len(seq) + 1, dtype=int)
    else:
        target_pos = np.asarray(parse_range(oligo_range))
    contigs_start = np.hstack([[target_pos[0]], target_pos[np.where(target_pos[1:] - target_pos[:-1] > 1)[0] + 1]])
    contigs_end = np.hstack([target_pos[np.where(target_pos[1:] - target_pos[:-1] > 1)], [target_pos[-1]]])
    assert len(contigs_start) == len(contigs_end)
    tile_sets = list()
    primer_pairs = list()
    for contig_start, contig_end in zip(contigs_start, contigs_end):
        primer_span_max_trial = primer_span_max
        complete_primer_set = False
        while primer_span_max_trial >= tile_len_min and not complete_primer_set:
            primer_pairs_trial = list()
            tile_sets_trial = list()
            if contig_start == 1:
                current_position = contig_start
            else:
                current_position = contig_start - (gap + 1)
            while True:
                # create forward primer
                direction = "dn" if current_position == 1 else "up"
                fwd_primer = grow_primer(seq, current_position, "fwd", direction, primer_tm_min)
                # create reverse primer
                if direction == "up":
                    fwd_primer_up, fwd_primer_dn = current_position - len(fwd_primer) + 1, current_position
                    current_position = min(current_position + primer_span_max_trial - len(fwd_primer), contig_end)
                elif direction == "dn":
                    fwd_primer_up, fwd_primer_dn = current_position, current_position + len(fwd_primer) - 1
                    current_position = min(current_position + primer_span_max_trial - 1, contig_end)
                else:
                    raise ValueError(direction)
                rev_primer = grow_primer(seq, current_position, "rev", "up", primer_tm_min)
                rev_primer_up, rev_primer_dn = current_position - len(rev_primer) + 1, current_position
                primer_pairs_trial.append((fwd_primer, rev_primer))
                # create tiles
                if fwd_primer_up == contig_start:
                    tile_up = fwd_primer_up
                else:
                    tile_up = fwd_primer_dn + (gap + 1)
                if rev_primer_dn == contig_end:
                    tile_dn = contig_end
                else:
                    tile_dn = rev_primer_up - (gap + 1)
                tile_span = tile_dn - tile_up + 1
                if tile_span < tile_len_min:
                    primer_span_max_trial -= 1
                    break
                else:
                    n_tiles = int(np.ceil(tile_span / tile_len_max))
                    tile_len_avg = tile_span / n_tiles
                    tile_ups = [tile_up + int(tile_len_avg * i) for i in range(n_tiles)]
                    tile_dns = [tile_ups[i] - 1 for i in range(1, n_tiles)] + [tile_dn]
                    tile_seqs = [revcomp(seq[up - 1: dn]) for up, dn in zip(tile_ups, tile_dns)]
                    tile_sets_trial.append(tile_seqs)
                if current_position == contig_end:
                    complete_primer_set = True
                    break
                current_position -= len(rev_primer) + gap * 2
        if complete_primer_set:
            primer_pairs.extend(primer_pairs_trial)
            tile_sets.extend(tile_sets_trial)
        else:
            raise ValueError("Failed")
    return primer_pairs, tile_sets


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("tile_len_min", type=int)
    parser.add_argument("tile_len_max", type=int)
    parser.add_argument("seq_len_max", type=int)
    parser.add_argument("primer_tm_min", type=float)
    parser.add_argument("--range")
    args = parser.parse_args()
    primer_pairs, tile_sets = olaygo(args.input, args.output, args.tile_len_min, args.tile_len_max, args.seq_len_max, args.primer_tm_min, args.range)
    print("\n".join([f"OC43-BAC_primer_{i}{d}\t{seq}\t" for i, seqs in enumerate(primer_pairs, start=1) for seq, d in zip(seqs, "FR")]))
    tiles = [tile for tile_set in tile_sets for tile in tile_set]
    print("\n".join([f"OC43-BAC_ASO_{i}-{j}\t{seq}\t" for i, seqs in enumerate(tile_sets, start=1) for j, seq in enumerate(seqs, start=1)]))
    print("Number of ASOs:", len([x for s in tile_sets for x in s]), "; lengths", min(map(len, tiles)), "-", max(map(len, tiles)))
    print("Number of primers:", len([x for s in primer_pairs for x in s]))

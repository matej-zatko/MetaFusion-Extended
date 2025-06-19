import csv
import os
import argparse
import sys

def parse_truth_set(truth_file):
    truth_fusions = set()
    with open(truth_file) as f:
        for line in f:
            sample, fusion = line.strip().split("|")
            genes = fusion.split("--")
            truth_fusions.add((sample, frozenset(genes)))
    return truth_fusions

def mark_tp_fp(input_file, truth_file, output_dir):
    try:
        os.makedirs(output_dir)
    except OSError:
        if not os.path.isdir(output_dir):
            raise

    output_file = os.path.join(output_dir, 'marked_fusions.tsv')
    fn_output = os.path.join(output_dir, 'false_negatives.txt')
    summary_output = os.path.join(output_dir, 'summary.txt')

    truth_fusions = parse_truth_set(truth_file)
    seen_tp = set()
    detected_fusions = set()

    counts = {'TP': 0, 'TP-DUP': 0, 'FP': 0}

    with open(input_file, 'rb') as infile, open(output_file, 'wb') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        fieldnames = reader.fieldnames + ['match']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for row in reader:
            sample = row['samples'].strip()
            genes1 = [g.strip() for g in row['gene1'].split(',') if g.strip()]
            genes2 = [g.strip() for g in row['gene2'].split(',') if g.strip()]

            match_label = 'FP'
            for g1 in genes1:
                for g2 in genes2:
                    if g1 == g2:
                        continue
                    fusion_key = (sample, frozenset([g1, g2]))
                    if fusion_key in truth_fusions:
                        detected_fusions.add(fusion_key)
                        if fusion_key not in seen_tp:
                            match_label = 'TP'
                            seen_tp.add(fusion_key)
                        else:
                            match_label = 'TP-DUP'
                        break
                if match_label != 'FP':
                    break

            counts[match_label] += 1
            row['match'] = match_label
            writer.writerow(row)

    # Write false negatives
    false_negatives = truth_fusions - detected_fusions
    with open(fn_output, 'w') as fn_file:
        for sample, genes in sorted(false_negatives):
            gene_list = sorted(genes)
            fn_file.write("%s|%s--%s\n" % (sample, gene_list[0], gene_list[1]))

    # Write summary statistics
    with open(summary_output, 'w') as sfile:
        sfile.write("TP: %d\n" % counts['TP'])
        sfile.write("TP-DUP: %d\n" % counts['TP-DUP'])
        sfile.write("FP: %d\n" % counts['FP'])
        sfile.write("FN: %d\n" % len(false_negatives))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mark fusions as TP, TP-DUP, FP; output FNs and summary.")
    parser.add_argument('--input', help='Input TSV file with fusion calls')
    parser.add_argument('--truth', help='Truth set file (format: sample|GENE1--GENE2)')
    parser.add_argument('--output_dir', help='Directory where outputs will be saved')
    args = parser.parse_args()

    if not args.input or not args.truth or not args.output_dir:
        parser.print_help()
        sys.exit(1)

    mark_tp_fp(args.input, args.truth, args.output_dir)
    print("Benchmarking results written to %s" % args.output_dir)
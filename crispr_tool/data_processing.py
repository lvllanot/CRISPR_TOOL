import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns



def reverse_complement(genome):
    complementary = {"A": "T", "C": "G", "G": "C", "T": "A"}
    r_complement = ""
    for base in genome[::-1]:
        r_complement += complementary[base]
    return r_complement


def process_fasta(filepath):

    entries = []
    current_header = ""
    current_sequence = ""

    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    parts = current_header.split()
                    header_first = parts[0]
                    code = header_first.split('.')[0]
                    name = " ".join(parts[1:])
                    entries.append({
                        "code": code,
                        "name": name,
                        "orientation": "forward",
                        "sequence": current_sequence.upper()
                    })
                current_header = line[1:]
                current_sequence = ""
            else:
                current_sequence += line

        if current_header:
            parts = current_header.split()
            header_first = parts[0]
            code = header_first.split('.')[0]
            name = " ".join(parts[1:]) if len(parts) > 1 else ""
            entries.append({
                "code": code,
                "name": name,
                "orientation": "forward",
                "sequence": current_sequence.upper()
            })

    df = pd.DataFrame(entries, columns=["code", "name", "orientation", "sequence"])
    df_reverse = df.copy()
    df_reverse['sequence'] = df_reverse['sequence'].apply(reverse_complement)
    df_reverse['orientation'] = 'reverse'
    df = pd.concat([df, df_reverse])
    return df


def run_local_blast(fasta_file, database, outname, blast_exe="blastn"):

    output_file = f"blast_out_{outname}.tsv"
    output_format = "6 qseqid pident qstart qend sstart send slen stitle"

    command = (f'{blast_exe} '
               f'-query {fasta_file} '
               f'-db {database} '
               f'-outfmt "{output_format}" '
               f'-out {output_file} '
               f'-task blastn-short '
               f'-word_size 20 '
               )
    os.system(command)
    return output_file


def database_to_fasta(df, output_filename="guides.fasta"):
    output_dir = os.path.dirname(output_filename)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_filename, "w") as f:
        for index, row in df.iterrows():
            header = f"> guide_{index}"
            sequence = row["sequence_guide"]
            f.write(header + "\n" + sequence + "\n")
    return output_filename


def plot_guides(total_length, df_guides, cas, section="all"):
    len_guide = 20
    sns.set(style="whitegrid", palette ="Greens_d")

    df_guides = df_guides.sort_values("position")

    fig, ax = plt.subplots(figsize=(12, 4))

    ax.hlines(y=0, xmin=0, xmax=total_length, linewidth=5, color = "#B0E0E6")
    ax.text(total_length / 2, -0.5, "Original Sequence", ha="center", va="top", fontsize=10, color="black")

    y_level = 2
    y_spacing = 2  # Espacio entre cada guía

    for idx, row in df_guides.iterrows():
        start = row["position"]
        end = start + len_guide

        guide_name = row.get("guide_name")

        ax.hlines(y=y_level, xmin=start, xmax=end, linewidth=4)
        ax.text((start + end) / 2, y_level, guide_name, ha="center", va="center", fontsize=5, color="black")

        y_level += y_spacing

    ax.set_xlim(0, total_length)
    ax.set_ylim(-1, y_level + 1)
    ax.set_xlabel("Position in sequence")
    ax.set_title(f"Guides of cas{cas} on Sequence ({section})")

    return


# total_length= 1640
# df_guides_example = "guides_cas9_all.tsv"
# cas = 9
# section = "all"
#
# plot_guides_lines(total_length, df_guides_example, cas, section)

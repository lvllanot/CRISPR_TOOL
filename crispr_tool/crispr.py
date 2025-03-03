import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

import data_processing as dp


def calculate_gc_content(seq):
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return round((gc_count / len(seq)) * 100, 2)


def design_guides(fasta_file, cas):

    output_dir = "../data"
    os.makedirs(output_dir, exist_ok=True)

    output_dir = "../results"
    os.makedirs(output_dir, exist_ok=True)

    df = dp.process_fasta(fasta_file)
    guides_list = []
    cas = str(cas)

    for index, row in df.iterrows():
        orientation = row["orientation"]
        seq = row["sequence"].upper()
        L = len(seq)

        if orientation == 'forward':
            if cas == "9":
                for i in range(L - 22):
                    candidate = seq[i:i + 20]
                    pam = seq[i + 20:i + 23]
                    GC = calculate_gc_content(candidate)
                    if len(pam) == 3 and pam[1:] == "GG":
                        guides_list.append({
                            "orientation": "forward",
                            "sequence_guide": candidate,
                            "position": i,
                            "GC" : GC
                        })
            elif cas == "12":
                for i in range(4, L - 20):
                    pam = seq[i - 4:i]
                    candidate = seq[i:i + 20]
                    GC = calculate_gc_content(candidate)
                    if pam[:3] == "TTT" and pam[3] in "ACG":
                        guides_list.append({
                            "orientation": "forward",
                            "sequence_guide": candidate,
                            "position": i,
                            "GC" :GC
                        })
        elif orientation == 'reverse':
            if cas == "9":
                for j in range(L - 22):
                    candidate = seq[j:j + 20]
                    pam = seq[j + 20:j + 23]
                    GC = calculate_gc_content(candidate)
                    if len(pam) == 3 and pam[1:] == "GG":
                        pos = L - j - 20
                        guides_list.append({
                            "orientation": "reverse",
                            "sequence_guide": candidate,
                            "position": pos,
                            "GC": GC
                        })
            elif cas == "12":
                for j in range(4, L - 20):
                    pam = seq[j - 4:j]
                    candidate = seq[j:j + 20]
                    GC = calculate_gc_content(candidate)
                    if pam[:3] == "TTT" and pam[3] in "ACG":
                        pos = L - j - 20
                        guides_list.append({
                            "orientation": "reverse",
                            "sequence_guide": candidate,
                            "position": pos,
                            "GC":GC})

    guides_df = pd.DataFrame(guides_list, columns=["orientation", "sequence_guide", "position", "GC"])

    guides_fasta = dp.database_to_fasta(guides_df, output_filename="guides.fasta")

    # verify_guides
    output_file_BD = f"../data/BD_editar"
    command_ver = (f'makeblastdb '
               f'-in {fasta_file} '
               f'-dbtype nucl '
               f'-out {output_file_BD} '
               )
    os.system(command_ver)

    blast_output_verify = dp.run_local_blast(guides_fasta, "../data/BD_editar", outname=f"guiasverificar")

    verify = pd.read_csv("blast_out_guiasverificar.tsv", sep="\t", names=['qseqid', 'pident', 'qstart', 'qend', 'sstart', 'send', 'slen', 'stitle'])
    verify = verify[['qseqid']]

    counts = verify['qseqid'].value_counts().reset_index()
    counts.columns = ['qseqid', 'hit_count']
    count_list = counts['hit_count'].tolist()
    name_list = counts['qseqid'].tolist()

    guides_df["guide_name"] = name_list
    guides_df["counts_gen"] = count_list
    guides_df = guides_df.iloc[:,[4, 0, 1, 2, 3, 5]]
    os.remove("blast_out_guiasverificar.tsv")
    # guides.to_csv("df_raw_guides.tsv", sep="\t")

    selection = []
    for index, row in guides_df.iterrows():
        if row["counts_gen"] == 1 and 40 <= row["GC"] <= 60:
            selection.append(1)
        else:
            selection.append(0)

    guides_df["selection"] = selection
    guides_selec = guides_df[guides_df["selection"] == 1]
    guides_selec = guides_selec.iloc[:, [0, 1, 2, 3, 4]]
    # guides_selec.to_csv("df_guides.tsv", sep="\t")

    return guides_selec


def crisprcas(fasta_file, cas, section = "all", graph = False ):
    guides_df = design_guides(fasta_file, cas)

    sequence = str()
    with open(fasta_file) as file:
        for line in file:
            if line.startswith(">"):
                pass
            else:
                sequence += line.strip()

    size = int(len(sequence))
    parts = int(size/3)
    selection = []

    if section == "all":
        pass
    elif section == "start":
        part_gen_s = 0
        part_gen_e = parts
        for index, row in guides_df.iterrows():
            if part_gen_s <= row["position"] < part_gen_e:
                selection.append(1)
            else:
                selection.append(0)
    elif section == "middle":
        part_gen_s = parts
        part_gen_e = 2 * parts
        for index, row in guides_df.iterrows():
            if part_gen_s <= row["position"] < part_gen_e:
                selection.append(1)
            else:
                selection.append(0)
    elif section == "end":
        part_gen_s = 2 * parts
        part_gen_e = size
        for index, row in guides_df.iterrows():
            if part_gen_s <= row["position"] <= part_gen_e:
                selection.append(1)
            else:
                selection.append(0)

    if section != "all":
        guides_df["selection"] = selection
        guides_df = guides_df[guides_df["selection"] == 1]

    guides_df = guides_df.iloc[:, [0, 1, 2, 3, 4]]
    out_name = f"../results/guides_cas{cas}_{section}.tsv"
    guides_df.to_csv(out_name, sep="\t", index=False)
    os.remove("guides.fasta")

    if graph:
        dp.plot_guides(size, guides_df, cas, section)
        plt.savefig(f"../results/guides_cas{cas}_{section}.svg")
    else:
        pass

    return


def crisprcas_next_to(fasta_file, cas, next_to = 0 , graph = False ):

    guides_df = design_guides(fasta_file, cas)

    sequence = str()
    with open(fasta_file) as file:
        for line in file:
            if line.startswith(">"):
                pass
            else:
                sequence += line.strip()

    size = int(len(sequence))

    if next_to < 0 or next_to > size:
        print("La posicion no esta dentro de la secuencia")
        sys.exit()
    else:
        pass


    selection = []

    part_gen_s = next_to - 100
    part_gen_e = next_to + 100
    for index, row in guides_df.iterrows():
        if part_gen_s <= row["position"] < part_gen_e:
            selection.append(1)
        else:
            selection.append(0)

    guides_df["selection"] = selection
    guides_df = guides_df[guides_df["selection"] == 1]

    guides_df = guides_df.iloc[:, [0, 1, 2, 3, 4]]
    out_name = f"../results/guides_cas{cas}_nextto_{next_to}.tsv"
    guides_df.to_csv(out_name, sep="\t", index=False)
    os.remove("guides.fasta")

    if graph:
        dp.plot_guides(size, guides_df, cas, next_to)
        plt.savefig(f"../results/guides_cas{cas}_nextto_{next_to}.svg")
        plt.show()
    else:
        pass

    return


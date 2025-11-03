# table_exporter.py
import csv

def save_frequency_table(freq_dicts: dict, output_file: str):
    genomes = list(freq_dicts.keys())
    keys = sorted(next(iter(freq_dicts.values())).keys())

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Genome"] + keys)
        for genome in genomes:
            row = [genome] + [f"{freq_dicts[genome][k]:.6f}" for k in keys]
            writer.writerow(row)
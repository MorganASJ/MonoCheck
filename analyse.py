import csv
from collections import defaultdict

def analyze_monophyly(csv_file):
    clade_status = defaultdict(lambda: {"monophyletic": 0, "total": 0})

    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            clade = row['Clade']
            status = row['Status']
            clade_status[clade]["total"] += 1
            if status == "monophyletic":
                clade_status[clade]["monophyletic"] += 1

    for clade, counts in clade_status.items():
        percentage = (counts["monophyletic"] / counts["total"]) * 100
        print(f"Clade: {clade}, Monophyletic: {counts['monophyletic']}, Total: {counts['total']}, Percentage: {percentage:.2f}%")

if __name__ == "__main__":
    csv_file = "tree_stats.csv"
    analyze_monophyly(csv_file)
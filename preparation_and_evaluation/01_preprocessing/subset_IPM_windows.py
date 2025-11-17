import pandas as pd

# Define the file path
file_path = "./../studies/deepCIS/ipmArthDAPmotifs/extracted_ranges_1.txt"
common_ids = pd.read_csv(file_path, sep=" ", header=None, names=["seq_id", "label"])
common_ids = common_ids.drop_duplicates(subset=["seq_id", "label"], keep="first")
selection = [
    "TCP",
    "WRKY",
    "bHLH",
    "BES1",
    "bZIP",
    "Trikelix",
    "MYB",
    "CAMPTRA",
    "FAR1",
    "MYB"
]


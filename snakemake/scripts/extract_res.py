import os
import re
from scipy.stats import chi2
import numpy as np

# Directory to scan
res_dir = "res/"

# Output table structure
header = ["NAME", "BRANCH", "REPL", "NULL_LOG-L", "ALT_LOG-L", "LRT", "p-value", "-logpvalue"]
data = {}

# Regex to extract relevant information from .res files
best_logl_pattern = re.compile(r"BEST LOG-L: ([\d\.\-]+)")

# Process each .res file
if os.path.isdir(res_dir):
    for filename in os.listdir(res_dir):
        if filename.endswith(".res"):
            filepath = os.path.join(res_dir, filename)
            with open(filepath, "r") as f:
                content = f.read()
                
                best_logl_match = best_logl_pattern.search(content)
                
                if best_logl_match:
                    best_logl = best_logl_match.group(1)
                    name_parts = filename.replace(".res", "").split(".")
                    name = ".".join(name_parts[:-3])
                    branch = name_parts[-3]
                    model = name_parts[-2]
                    replicate = name_parts[-1]
                    
                    key = (name, branch, replicate)
                    
                    if key not in data:
                        data[key] = {"null_logl": "NA", "alt_logl": "NA"}
                    
                    if model == "null":
                        data[key]["null_logl"] = best_logl
                    elif model == "alt":
                        data[key]["alt_logl"] = best_logl

    # Compute LRT, p-value, and -log(p-value)
    for key, values in data.items():
        if values["null_logl"] != "NA" and values["alt_logl"] != "NA":
            null_logl = float(values["null_logl"])
            alt_logl = float(values["alt_logl"])
            lrt = -2 * (null_logl - alt_logl)
            pval = 1 - chi2.cdf(lrt, 1)
            log_pval = -np.log(pval)
            values["lrt"] = lrt
            values["pval"] = pval
            values["log_pval"] = log_pval
        else:
            values["lrt"] = "NA"
            values["pval"] = "NA"
            values["log_pval"] = "NA"

    # Write the summary table
    with open("summary_table.txt", "w") as out_file:
        out_file.write("\t".join(header) + "\n")
        for (name, branch, replicate), values in sorted(data.items()):
            row = [name, branch, replicate, values["null_logl"], values["alt_logl"], str(values["lrt"]), str(values["pval"]), str(values["log_pval"])]
            out_file.write("\t".join(row) + "\n")

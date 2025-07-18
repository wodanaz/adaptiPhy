import re  # for regular expressions
import sys

# Get the input file (a file containing one or more file paths)
infile = sys.argv[1]

with open(infile) as f:
    reflist = f.read().splitlines()

for filepath in reflist:
    with open(filepath, 'r') as myfile:
        data = myfile.read().splitlines()
    
    # Open the output file with .prunned suffix
    output_filename = f"{filepath}.prunned"
    with open(output_filename, 'w') as output:
        results = {}
        seq = ""
        compName = ""
        
        # Parse the FASTA file
        for line in data:
            if line.startswith(">"):
                if compName != "":
                    results[compName] = seq
                compName = line[1:]
                seq = ""
            else:
                seq = seq + line
        results[compName] = seq
        
        # Calculate sequence length using the first sequence
        nsize = len(results[list(results.keys())[0]])
        # Gather positions where the base is not one of A, C, G, T
        indels = set()
        for spec in results:
            for i in range(nsize):
                base = results[spec][i]
                if base.upper() not in ["A", "C", "G", "T"]:
                    indels.add(i)
        
        # Create prunned sequences by removing positions with indels
        prunned = {}
        for spec in results:
            new_seq = ""
            for i in range(nsize):
                if i not in indels:
                    new_seq = new_seq + results[spec][i]
            prunned[spec] = new_seq
        
        # Write the prunned sequences to the output file
        for j, spec in enumerate(list(prunned.keys())):
            if j > 0:
                output.write("\n")
            output.write(">" + spec + "\n")
            output.write(prunned[spec])

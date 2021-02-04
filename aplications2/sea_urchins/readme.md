# Workflow for analysis of evolution by positive selection in sea urchins


Measure local compositional complexity (LCC) of DNA sequence


```bash 

nano dolcc_content.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
python LCC_content.py


sbatch dolcc_content.sh

```

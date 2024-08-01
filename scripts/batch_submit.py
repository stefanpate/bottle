import subprocess
from dataclasses import dataclass

@dataclass
class BatchScript:
    '''
    allocation:str
        p30041 | b1039
    partition:str
        gengpu | short | normal | long | b1039
    mem:str
        memory required e.g., 8G
    time:str
        hours of compute e.g., 4
    file:str
        script name, e.g., file.py
    arg_str:
        string of args following script name e.g,. -d 4 --example
    '''
    allocation:str
    partition:str
    mem:str
    n_cores:int
    time:str
    script:str

    def write(self, arg_str, job_name):
        cpu_blacklist = ["#SBATCH --gres=gpu:a100:1"]
        lines = [
            f"#!/bin/bash",
            f"#SBATCH -A {self.allocation}",
            f"#SBATCH -p {self.partition}",
            f"#SBATCH --gres=gpu:a100:1",
            f"#SBATCH -N 1",
            f"#SBATCH -n {self.n_cores}",
            f"#SBATCH --mem={self.mem}",
            f"#SBATCH -t {self.time}:00:00",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --output=../logs/out/{job_name}",
            f"#SBATCH --error=../logs/error/{job_name}",
            f"#SBATCH --mail-type=END",
            f"#SBATCH --mail-type=FAIL",
            f"#SBATCH --mail-user=stefan.pate@northwestern.edu",
            f"ulimit -c 0",
            f"module load python/anaconda3.6",
            f"source activate mine",
            f"python -u {self.script} {arg_str}",
        ]

        if self.partition == 'gengpu':
            return '\n'.join(lines)
        else:
            return '\n'.join([line for line in lines if line not in cpu_blacklist])


if __name__ == '__main__':
    from src.utils import load_json, save_json
    allocation = 'b1039'
    partition = 'b1039'
    mem = '10GB'
    n_cores = 1
    hours = 12
    n_chunks = 20
    full = load_json("/home/spn1560/bottle/data/sprhea/sprhea_240310_v3.json")
    batch_script = BatchScript(
        allocation=allocation,
        partition=partition,
        mem=mem,
        n_cores=n_cores,
        time=hours,
        script='map_operators.py'
    )
    
    split_up = [{} for _ in range(n_chunks)]
    for idx, key in enumerate(full):
        sub_idx = idx % n_chunks

        split_up[sub_idx][key] = full[key]

    for i in range(n_chunks):
        save_json(split_up[i], f"../data/sprhea/sprhea_v3_part_{i}.json")
        

    for i in range(n_chunks):
        arg_str = f"../data/rules/minimal1224_all_uniprot.tsv ../data/sprhea/sprhea_v3_part_{i}.json ../artifacts/operator_mapping_sprhea_v3_min_ops_part_{i}.tsv"
        job_name = f"sprhmop{i}"
        shell_script = batch_script.write(arg_str, job_name)
                    
        with open("batch.sh", 'w') as f:
            f.write(shell_script)

        subprocess.run(["sbatch", "batch.sh"])

'''
e.g.
python -u map_operators.py /home/spn1560/bottle/data/rules/minimal1224_all_uniprot.tsv 
/home/spn1560/bottle/data/sprhea/sprhea_v3_part_19.json 
/home/spn1560/bottle/artifacts/operator_mapping_sprhea_v3_min_ops_part_19.tsv
'''

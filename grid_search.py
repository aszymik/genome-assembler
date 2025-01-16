import wandb
import assembly
import subprocess

# Weak:
# reads_coverage
# fragmentation_score

def evaluate_model(config=None):
    """
    Evaluate the model with the given parameters.
    
    :param config: A dictionary containing the parameters for the run.
    """
    with wandb.init(config=config):
        config = wandb.config

        # Parameters from the sweep
        k = config.k
        kmer_thresh = config.kmer_thresh
        weight_thresh = config.weight_thresh
        # similarity_thresh = config.similarity_thresh
        similarity_thresh = 3

        
        # Evaluate the model
        scores = get_scores(k, kmer_thresh, weight_thresh, similarity_thresh)

        # Log metrics
        wandb.log({
            'ref_coverage': scores[0], 
            'reads_coverage': scores[1], 
            'identity_score': scores[2], 
            'fragmentation_score': scores[3],
            'score': scores[4]
        })


def get_scores(k, kmer_thresh, weight_thresh, similarity_thresh):
    """
    Run the assembly program and calculate score.
    """
    assembly.main('training/reads/reads1.fasta', 'outs/reads1_contigs', k=k, kmer_thresh=kmer_thresh, weight_thresh=weight_thresh, similarity_thresh=similarity_thresh)
    result = subprocess.run(['./evaluate.sh', 'outs/reads1_contigs'], capture_output=True)

    # Decode subprocess output
    output = result.stdout.decode('utf-8')
    lines = output.strip().split('\n')

    scores = []
    for line in lines:
        if ':' in line:
            _, score = line.split(':', 1)
            scores.append(score.strip())
    return scores



sweep_config = {
    "method": "bayes",
    "metric": {"name": "score", "goal": "maximize"},
    "parameters": {
        "k": {"values": list(range(13, 37, 2))},
        "kmer_thresh": {"values": list(range(1, 8))},
        "weight_thresh": {"values": list(range(1, 6))},
        # "similarity_thresh": {"values": list(range(1, 11))},
    }
}

wandb.login()
sweep_id = wandb.sweep(sweep_config, project='genome-assembler')
wandb.agent(sweep_id, function=evaluate_model)



# run_gene_flux.py
import os, argparse, json
import numpy as np
import matplotlib.pyplot as plt
from gene_flux_engine import run_geneflux
from sequence_tools import codon_counts, aa_counts

parser = argparse.ArgumentParser()
parser.add_argument("--generations", type=int, default=50)
parser.add_argument("--pop", type=int, default=40)
parser.add_argument("--mutprob", type=float, default=0.25)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

OUT_GRAPHS = "../graphs"
OUT_LOGS = "../logs"
os.makedirs(OUT_GRAPHS, exist_ok=True)
os.makedirs(OUT_LOGS, exist_ok=True)

# default initial sequence (short example); replace with your FASTA if desired
initial = ["ATGGCTGCTGCTGAAATGGCATGCTGCTAGCTGACTATGCATGCGT"]
res = run_geneflux(initial_seqs=initial,
                   generations=args.generations,
                   pop_size=args.pop,
                   mutation_prob=args.mutprob,
                   seed=args.seed,
                   out_dir=OUT_GRAPHS,
                   log_path=os.path.join(OUT_LOGS, "geneflux_log.json"))

history = res["history"]
gens = [h["generation"] for h in history]
avg = [h["avg_score"] for h in history]
mx = [h["max_score"] for h in history]

# fitness plot
plt.figure(figsize=(7,4))
plt.plot(gens, avg, label="avg_score")
plt.plot(gens, mx, label="max_score")
plt.xlabel("Generation")
plt.ylabel("Score")
plt.title("GeneFlux: Fitness over Generations")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_GRAPHS, "geneflux_fitness.png"))
plt.close()

# build aggregated codon frequency heatmap (generations x codon positions)
# choose a fixed codon order
codon_order = sorted([a+b+c for a in "ATGC" for b in "ATGC" for c in "ATGC"])
max_codons = 0
for h in history:
    cc = h.get("codon_counts_best", {})
    max_codons = max(max_codons, sum(cc.values()))

heat = np.zeros((len(history), len(codon_order)))
for i,h in enumerate(history):
    cc = h.get("codon_counts_best", {})
    total = sum(cc.values()) or 1
    for j,c in enumerate(codon_order):
        heat[i,j] = cc.get(c,0)/total

plt.figure(figsize=(10,5))
plt.imshow(heat.T, aspect='auto', origin='lower', interpolation='nearest')
plt.xlabel("Generation")
plt.ylabel("Codon (sorted)")
plt.title("GeneFlux: Codon frequency (best sequence per generation)")
plt.colorbar(label="frequency")
plt.tight_layout()
plt.savefig(os.path.join(OUT_GRAPHS, "geneflux_codon_heatmap.png"))
plt.close()

# final best sequence save
final_best = res["final_best"]
with open(os.path.join(OUT_GRAPHS, "geneflux_final_best_sequence.txt"), "w") as f:
    f.write(final_best)

print("Done. Graphs saved to:", os.path.abspath(OUT_GRAPHS))
print("Log saved to", os.path.abspath(os.path.join(OUT_LOGS, "geneflux_log.json")))

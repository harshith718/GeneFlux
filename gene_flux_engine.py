# gene_flux_engine.py
# Core engine: population of sequences, mutation operators, selection, and logging

import random, statistics, json, os
from typing import List, Dict, Any, Tuple
from sequence_tools import clean_seq, translate, codon_counts, gc_content

BASES = ["A","T","G","C"]

def random_substitution(seq: str) -> Tuple[str,int,str]:
    if not seq:
        return seq, -1, ''
    pos = random.randrange(len(seq))
    old = seq[pos]
    choices = [b for b in BASES if b != old]
    new = random.choice(choices)
    new_seq = seq[:pos] + new + seq[pos+1:]
    return new_seq, pos, new

def random_insertion(seq: str) -> Tuple[str,int,str]:
    pos = random.randrange(len(seq)+1)
    b = random.choice(BASES)
    return seq[:pos] + b + seq[pos:], pos, b

def random_deletion(seq: str) -> Tuple[str,int,str]:
    if len(seq) <= 1:
        return seq, -1, ''
    pos = random.randrange(len(seq))
    return seq[:pos] + seq[pos+1:], pos, seq[pos]

def apply_mutation(seq: str, mut_type: str) -> Dict[str,Any]:
    if mut_type == 'sub':
        s,pos,b = random_substitution(seq)
        return {"seq": s, "type": "sub", "pos": pos, "base": b}
    if mut_type == 'ins':
        s,pos,b = random_insertion(seq)
        return {"seq": s, "type": "ins", "pos": pos, "base": b}
    if mut_type == 'del':
        s,pos,b = random_deletion(seq)
        return {"seq": s, "type": "del", "pos": pos, "base": b}
    return {"seq": seq, "type":"none", "pos": -1, "base": ""}

def mutate(seq: str, mutation_prob: float, mutation_rates: Dict[str,float]) -> Dict[str,Any]:
    if random.random() > mutation_prob:
        return {"seq": seq, "type":"none", "pos":-1, "base":""}
    r = random.random()
    cum = 0.0
    for mt,rate in mutation_rates.items():
        cum += rate
        if r <= cum:
            return apply_mutation(seq, mt)
    return apply_mutation(seq, "sub")

def score_sequence(seq: str) -> float:
    # Composite score using GC closeness to 50% and translated stop-codon penalty
    s = clean_seq(seq)
    if not s:
        return 0.0
    gc = gc_content(s)/100.0
    gc_score = 1.0 - abs(gc - 0.5)  # 1.0 when GC == 0.5
    aa = translate(s)
    stops = aa.count('*')
    stop_penalty = max(0.0, 1.0 - 0.5 * (stops))
    length_factor = max(0.2, 1.0 - abs(len(s) - 120)/500.0)
    return max(0.0, gc_score * stop_penalty * length_factor)

class Population:
    def __init__(self, sequences: List[str]):
        self.seqs = [clean_seq(s) for s in sequences]

    def evaluate(self) -> List[float]:
        return [score_sequence(s) for s in self.seqs]

    def select(self, frac: float) -> List[str]:
        scored = list(zip(self.seqs, self.evaluate()))
        scored.sort(key=lambda x: x[1], reverse=True)
        keep = max(1, int(len(scored) * frac))
        return [s for s,_ in scored[:keep]]

def run_geneflux(initial_seqs: List[str],
                 generations:int=50,
                 pop_size:int=40,
                 mutation_prob:float=0.2,
                 mutation_rates:Dict[str,float]=None,
                 selection_frac:float=0.5,
                 seed:int=42,
                 out_dir:str="../graphs",
                 log_path:str="../logs/geneflux_log.json") -> Dict[str,Any]:
    if mutation_rates is None:
        mutation_rates = {"sub":0.7, "ins":0.15, "del":0.15}
    random.seed(seed)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    pop = (initial_seqs * ((pop_size//len(initial_seqs))+1))[:pop_size]
    population = Population(pop)
    history = []

    for g in range(generations):
        new_seqs = []
        gen_mutations = []
        for seq in population.seqs:
            res = mutate(seq, mutation_prob, mutation_rates)
            new_seqs.append(res["seq"])
            gen_mutations.append(res)
        population.seqs = new_seqs

        scores = population.evaluate()
        avg = sum(scores)/len(scores)
        mx = max(scores)
        best = population.seqs[scores.index(mx)]

        # selection & refill
        survivors = population.select(selection_frac)
        new_pop = survivors.copy()
        while len(new_pop) < pop_size:
            parent = random.choice(survivors)
            new_pop.append(parent)
        population.seqs = new_pop

        # record codon / aa freq for best
        codon_map = codon_counts(best)
        aa = translate(best)
        history.append({
            "generation": g,
            "avg_score": avg,
            "max_score": mx,
            "best_seq": best,
            "best_aa": aa,
            "mutations_sample": gen_mutations[:8],
            "codon_counts_best": codon_map
        })

    with open(log_path, "w") as f:
        json.dump({"history": history}, f, indent=2)

    return {"generations": generations, "history": history, "final_best": history[-1]["best_seq"]}

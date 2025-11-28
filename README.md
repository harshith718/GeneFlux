# GeneFlux — Genetic Drift & Adaptive Dynamics Simulator

GeneFlux is a lightweight evolutionary modelling project focused on understanding **genetic drift**, **allele-frequency changes**, and **adaptive shifts** under different evolutionary pressures.  
It simulates neutral evolution, selection-driven evolution, and hybrid models that combine mutation, drift, and fitness constraints.

GeneFlux is part of a 6-project computational evolution research suite exploring mutation, selection, speciation, enzyme design, and ecological adaptation.

---

## 🔬 Overview

GeneFlux answers fundamental evolutionary questions:

- How do allele frequencies shift over generations?
- What happens when drift dominates over selection?
- How do beneficial mutations compete with random fluctuations?
- How does population size influence evolutionary noise?
- What trajectories emerge under neutral vs selective simulations?

It provides a clear, visual, and code-driven way to explore these patterns.

---

## 📁 Project Structure

```
GeneFlux/
│
├── code/
│   ├── geneflux_simulator.py
│   ├── allele_drift_engine.py
│   ├── selection_model.py
│   ├── run_geneflux_plot.py
│   └── example_allele_dataset.json
│
├── graphs/
│   ├── allele_frequency_curve.png
│   ├── drift_vs_selection.png
│   └── population_noise_plot.png
│
└── logs/
    ├── drift_run_log.json
    └── best_allele_trajectory.txt
```

- **code/** → All Python modules and engines  
- **graphs/** → Visualization outputs  
- **logs/** → Logs, run outputs, and final state summaries  

---

## ▶️ How to Run

### **1. Install Python (3.9+)**
Check your version:
```bash
python --version
```

### **2. Install dependencies**
```
pip install numpy matplotlib
```

### **3. Run a GeneFlux simulation**
```
python code/run_geneflux_plot.py
```

This will:

- simulate drift or drift+selection  
- track allele pools  
- generate graphs  
- save the best trajectory  
- write logs for reproducibility  

Outputs appear in:

- `graphs/`
- `logs/`

---

## 📊 Outputs Generated

- **allele_frequency_curve.png** — dominant allele frequency across generations  
- **drift_vs_selection.png** — comparison of evolutionary noise vs selection effects  
- **population_noise_plot.png** — visualizing random drift intensity  
- **best_allele_trajectory.txt** — most successful allele trajectory  
- **drift_run_log.json** — full simulation log including parameters and steps  

---

## 🎯 Purpose

GeneFlux is designed to:

- teach core population genetics concepts  
- provide computational intuition of drift vs selection  
- demonstrate simulation-based evolutionary analysis  
- support deeper portfolio work in evolutionary computation  

It is optimized for clarity, reproducibility, and undergraduate-level research.

## 🔗 Portfolio Link  
Complete 6-project evolution research collection:  
https://west-route-a3b.notion.site/BioGraph-Evolution-Research-Portfolio-2b69325d1ab1804dab15f731b8af6581?source=copy_link

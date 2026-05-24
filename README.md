# Replication material for the ESE Heckman selection paper

This repository contains the computational material for the paper:

**Transformation--Conditioning Heckman-Type Selection Models with Extended Skew-Elliptical Errors**  
Roberto Vila, Cecilia Castro, Víctor Leiva, Helton Saulo, Carolina Marchant.

The repository provides R code and selected numerical outputs for the extended skew-elliptical (ESE) framework for Heckman-type sample-selection models. The paper develops a transformation--conditioning, or hidden-truncation, construction for bivariate sample-selection models and studies its likelihood-based implementation through Monte Carlo recovery experiments, synthetic benchmark comparisons, and an empirical application.

The computational material covers:

- correctly specified ESN/EST Monte Carlo recovery studies;
- synthetic benchmark comparisons against normal Heckman, Student--\(t\) Heckman, and contaminated-normal Heckman alternatives;
- an empirical application to the classical `Mroz87` dataset;
- scenario-level summaries, replication-level outputs, tables, figures, metadata, random-seed settings, and R session information;
- scripts designed for local execution and high-performance execution through Slurm arrays.

The benchmark results should be interpreted as a likelihood--penalty comparison across regimes. The proposed ESN and EST specifications can improve the unpenalized likelihood in some hidden-truncation designs, but the improvement may be too small to offset the additional complexity under AIC or BIC in finite samples. The repository therefore supports the interpretation of the ESE Heckman class as an expanded modeling toolkit to be assessed through likelihood, information criteria, parsimony, convergence diagnostics, boundary diagnostics, and computational cost.

## Repository structure

```text
scripts/
  R scripts used for Monte Carlo recovery experiments, synthetic benchmark
  comparisons, empirical analysis, and report generation.

results/
  Generated outputs from the simulation and empirical workflows. These folders
  may include per-replication CSV files, scenario summaries, report tables,
  figures, metadata, R session information, and serialized R objects.

data/
  Input data or selected summary CSV files used for reproducing tables and
  figures in the paper.

logs/
  Slurm output and error logs for the high-performance runs.

downloads/
  Optional compressed archives prepared for transfer from the HPC system.
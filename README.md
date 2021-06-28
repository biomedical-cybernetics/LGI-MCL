# Latent Geometry Inspired - Markov Clustering (LGI-MCL)

Here you will find the scripts used to compute the results of the ANS article:

## Geometrical Inspired Pre-Weighting Enhances Markov Clustering Community Detection in Complex Networks

**Abstract** — Markov clustering is an effective unsupervised pattern recognition algorithm for data clustering in high-dimensional feature space. However, its community detection performance in complex networks has been demonstrating results far from the state of the art methods such as Infomap and Louvain. The crucial issue is to convert the unweighted network topology in a ‘smart-enough’ pre-weighted connectivity that adequately steers the stochastic flow procedure behind Markov clustering. Here we introduce a conceptual innovation and we discuss how to leverage network latent geometry notions in order to design similarity measures for pre-weighting the adjacency matrix used in Markov clustering community detection. Our results demonstrate that the proposed strategy improves Markov clustering significantly, to the extent that it is often close to the performance of current state of the art methods for community detection. These findings emerge considering both synthetic ‘realistic’ networks (with known ground-truth communities) and real networks (with community metadata), and even when the real network connectivity is corrupted by noise artificially induced by missing or spurious links. Our study enhances the generalized understanding of how network geometry plays a fundamental role in the design of algorithms based on network navigability.

Please cite: Durán, C., Muscoloni, A. & Cannistraci, C.V. Geometrical inspired pre-weighting enhances Markov clustering community detection in complex networks. Appl Netw Sci 6, 29 (2021). https://doi.org/10.1007/s41109-021-00370-x

## Simulations:

### nPSO 

This folder contains the scripts to generate the nPSO network: nPSOmdel_networks
and corresponding MCL, Infomap and Louvain scripts to run on the nPSO networks.

### real

This folder contains the real networks analyzed in the manuscript
and corresponding MCL, Infomap and Louvain scripts to run the simulations

### real_perturbed

This folder contains the perturbed real networks analyzed in the manuscript, both with additional and removed links
and corresponding MCL, Infomap and Louvain scripts to run the simulations


## Note

The simulations were applied in a Windows environment. Hence some adjustments might be required to run them such as:
Previous installation of MCL and Infomap.
Previous installation of R and respective packages to run Louvain.
Adjust paths inside the scripts run_commdet_[simulationToRun].m

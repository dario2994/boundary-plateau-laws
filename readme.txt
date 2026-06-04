Pipeline:

a) ./enumerate_twice_derived
It requires >= 16gb of ram, and approximately six hours to be executed.
It creates data/twice_derived_computer.txt and data/twice_derived_human.txt.
These are lists of candidates for being twice-derived nets.
Each element of the list is a spherical map.

b) ./enumerate_once_derived
It requires a few minutes to be executed.
It reads the list of twice-derived nets from data/twice_derived_computer.txt
and creates data/once_derived_computer.txt, data/once_derived_human.txt,
data/once_derived_for_visualization.txt.
These are lists of candidates for being once-derived nets.
Each element of the list is a marked spherical map.

c) ./embed_or_exclude
Requires 2 hours to be executed.
It reads the list of once-derived nets from data/once_derived_computer.txt.
For each of them, it tries to embed it in the sphere.
It either manages to embed it (up to error ~0.00001) or proves that it
cannot be embedded. It also proves (for most nets) the (approximate) rigidity
of the embedding.
The ids of the embeddable nets are saved in data/embeddable_nets.txt
For each embeddable net, the data of a valid embedding is saved in
data/metrization/net_id.txt.

e) python scripts/run_visualization net_id
It visualizes the marked combinatorial map (in the once-derived nets list)
with net_id as well as an embedding in the sphere (if one was found).

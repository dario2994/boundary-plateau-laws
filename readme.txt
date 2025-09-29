Pipeline:

a) ./enumerate_derived_twice
It requires >= 16gb of ram, and a few hours to be executed.
It creates data/derived_twice_computer.txt and data/derived_twice_human.txt.
These are lists of candidates for being derived-twice nets.

b) ./enumerate_derived_once
It requires a few minutes to be executed.
It reads the list of derived-twice nets from data/derived_twice_computer.txt
and creates data/derived_once_computer.txt, data/derived_once_human.txt,
data/derived_once_for_visualization.txt.
These are lists of candidates for being derived-once nets.

c) ./embed_or_exclude_nets
Requires ??? hours to be executed.
It reads the list of derived-once nets from data/derived_once_computer.txt.
For each of them, it tries to embed it in the sphere.
It either manages to embed it (up to error ~0.00001) or proves that it
cannot be embedded. It also proves (for most nets) the (approximate) rigidity
of the embedding.
The ids of the embeddable nets are saved in data/embeddable_nets.txt
For each embeddable net, the data of a valid embedding is saved in data/metrization/net_id.txt.

e) python scripts/run_visualization net_id
It visualizes the "combinatorial" derived once net with net_id as well
as an embedding in the sphere (if one was found).

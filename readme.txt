Pipeline:

a) ./main_combinatorial_classification
It requires >= 16gb of ram, and a few hours to be executed.
It creates data/derived_twice_computer.txt and data/derived_twice_human.txt.
These are lists of candidates for being derived-twice nets.

b) ./main_geometric_classification
It requires a few minutes to be executed.
It reads the list of derived-twice nets from data/derived_twice_computer.txt
and creates data/derived_once_computer.txt, data/derived_once_human.txt,
data/derived_once_for_visualization.txt.
These are lists of candidates for being once-twice nets.
         
c) ./sieve_derived_once_nets
It requires a fraction of a second to be executed.
It reads the list of derived-once nets from data/derived_once_computer.txt
and applies a sequence of criteria to prove that some nets are not realizable.
The ids of the nets that survive are saved in data/derived_once_survivors.txt

d) ./embed_or_exclude_nets
Requires ~60 hours to be executed.
It reads the list of survived derived-once nets from data/derived_once_survivors.txt.
For each of them, it tries to embed it in the sphere.
It either manages to embed it (up to error ~0.00001) or proves that it 
cannot be embedded.
There are a few nets (see the comments in the main) that need to be
excluded through other means.
The ids of the embeddable nets are saved in data/derived_once_metric_survivors.txt
For each embeddable net, the data of a valid embedding is saved in data/metrization/net_id.txt.

e) python run_visualization net_id
It visualizes the "combinatorial" derived once net with net_id as well
as an embedding in the sphere (if one was found).

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

c) python visualize_derived_once_nets.py
It requires ten minutes to be executed.
It reads the list of derived-once nets from data/derived_once_for_visualization.txt
and, for each derived-once net, creates data/derived_once_images/id.png
where id is the id of the net.
ACHTUNG: It reads from data/derived_once_for_visualization.txt as the
         python library we are using for planar embedding does not 
         support multiple edges and thus we tweak the nets so that
         they don't have multiple edges.
         
d) ./sieve_derived_once_nets
It requires a fraction of a second to be executed.
It reads the list of derived-once nets from data/derived_once_computer.txt
and applies a sequence of criteria to prove that some nets are not realizable.
The ids of the nets that survive are saved in data/derived_once_survivors.txt

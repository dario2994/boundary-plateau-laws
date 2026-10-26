Pipeline:

1) ./enumerate_twice_derived
It requires >= 16gb of ram, and approximately six hours to be executed.
It creates data/twice_derived_computer.txt and data/twice_derived_human.txt.
These are lists of candidates for being twice-derived nets.
Each element of the list is a spherical graph.

2) ./enumerate_once_derived
It requires a few minutes to be executed.
It reads the list of twice-derived nets from data/twice_derived_computer.txt
and creates data/once_derived_computer.txt, data/once_derived_human.txt,
data/once_derived_for_visualization.txt.
These are lists of candidates for being once-derived nets.
Each element of the list is a marked spherical graph.

3) ./embed_or_exclude
Requires 2 hours to be executed.
It reads the list of once-derived nets from data/once_derived_computer.txt.
For each of them, it tries to embed it in the sphere.
It either manages to embed it (with a given tolerance) or proves that it
cannot be embedded. It also proves (for most nets) the (approximate) rigidity
of the embedding.
The ids of the embeddable nets are saved in data/embeddable_nets.txt
For each embeddable net, the data of a valid embedding is saved in
data/metrized_nets/net_id.txt.

4) ./find_competitor
For each family of embeddable nets that must be excluded (rigid nets,
rigid nets with flexible pins, goggles' family, and fat fish family);
it finds a competitor (or many competitors for families with 1 or 2 parameters)
using Brakke's surface evolver and proves that no element of the family can be
a minimizer.
It saves the competitors for family id in data/competitors/id/.
For this to run, you must have Brakke's surface-evolver installed with the
patch:
https://gist.githubusercontent.com/dario2994/25b6b268ffd76eeef8e8d36904a12c09/raw/c47eb0c5fe9b479430bb17be794d1fb20a748c07/guide.md

Additional files:

A) python website/generate_site.py
Generates the data necessary to visualize the nets in the browser.
Should be run after the four main step have completed.
Then, run python3 -m http.server 8000 in the website directory and open the
browser at localhost:8000

B) [Not useful if you can do A)]
    python scripts/visualize_once_derived_topology.py --many --in data/once_derived_for_visualization.txt --out OUTPUT_DIRECTORY
Generates an .svg image for each once-derived topology (a marked combinatorial
map)  and saves it in the OUTPUT_DIRECTORY.

C) [Not useful if you can do A)]
    python scripts/visualize_once_derived_embedded.py --in data/metrized_nets/ID.txt --out OUTPUT_FILE.html
Generates a .html file (that you can then open with your browser) visualizing
the given once-derived net embedded in the sphere.

D) [Not useful if you can do A)]
    python scripts/visualize_competitor.py --in data/competitors/ID/COMPETITOR_FILE.txt --out OUTPUT_FILE.html
Generates an .html file (that you can then open with your browser) visualizing
the given competitor.

E) ./pin-lower-bound
Verifies the positivity of a function, used in the proof of
Proposition 11.2 (cf. Appendix F.1)

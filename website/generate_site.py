import pathlib
import json
import subprocess
import sys

site_dir = pathlib.Path(__file__).resolve().parent
data_dir = site_dir.parent / "data"
scripts_dir = site_dir.parent / "scripts"

if str(scripts_dir) not in sys.path:
    sys.path.append(str(scripts_dir))

from visualize_once_derived_embedded import parse_metrization_file, build_3d_net_traces
from visualize_competitor import parse_competitor_file


def get_net_statuses():
    nets = {}
    competitors_map = {}

    # 1. Read embeddable net IDs safely
    embeddable_ids = set()
    embeddable_file = data_dir / "embeddable_nets.txt"
    if embeddable_file.exists():
        with open(embeddable_file, "r") as f:
            content = f.read().split()
            if content:
                tokens = [t for t in content if t.isdigit()]
                if tokens and int(tokens[0]) == len(tokens) - 1:
                    tokens = tokens[1:]
                embeddable_ids = set(tokens)

    # 2. Get total once-derived count from once_derived_computer.txt
    total_nets = -1
    once_derived_file = data_dir / "once_derived_computer.txt"
    if once_derived_file.exists():
        with open(once_derived_file, "r") as f:
            first_line = f.readline().strip()
            if first_line.isdigit():
                total_nets = int(first_line)

    competitors_dir = data_dir / "competitors"
    special_red_ids = {158, 681, 920, 1100}

    # 3. Evaluate all net IDs explicitly
    for net_id in range(total_nets):
        has_metrization = str(net_id) in embeddable_ids

        net_comp_dir = competitors_dir / str(net_id)
        if net_comp_dir.exists() and net_comp_dir.is_dir():
            all_files = [p.name for p in net_comp_dir.glob(f"{net_id}*_competitor.txt")]

            valid_files = []
            for fname in all_files:
                clean = fname.replace("_competitor.txt", "")
                parts = clean.split("_")
                assert int(parts[0]) == net_id
                if net_id == 827 and len(parts) != 5:
                    continue
                try:
                    [float(p) for p in parts[1:]]
                    valid_files.append(fname)
                except ValueError:
                    continue

            if valid_files:
                param_count = (len(valid_files[0].split("_")) - 2) // 2
                competitors_map[net_id] = {
                    "files": sorted(valid_files),
                    "param_count": param_count,
                }

        # Status assignment precedence
        if net_id in special_red_ids:
            nets[net_id] = "red"
        elif has_metrization or (net_id in competitors_map):
            nets[net_id] = "green"
        else:
            nets[net_id] = "gray"

    return nets, competitors_map, embeddable_ids


# 1. Generate 2D topology SVGs in website/once_derived_topology/
topology_input = data_dir / "once_derived_for_visualization.txt"
topology_out_dir = site_dir / "once_derived_topology"
if topology_input.exists():
    subprocess.run(
        [
            sys.executable,
            str(scripts_dir / "visualize_once_derived_topology.py"),
            "--many",
            "--in",
            str(topology_input),
            "--out",
            str(topology_out_dir),
        ],
        check=True,
    )

nets, competitors_map, embeddable_ids = get_net_statuses()

# 2. Generate 3D JSON trace data for embeddable nets in website/once_derived_embedded/
embedded_out_dir = site_dir / "once_derived_embedded"
embedded_out_dir.mkdir(parents=True, exist_ok=True)

available_3d_ids = set()
for net_id_str in embeddable_ids:
    metrized_txt = data_dir / "metrized_nets" / f"{net_id_str}_visualization.txt"
    if metrized_txt.exists():
        try:
            verts, edges = parse_metrization_file(metrized_txt)
            traces = build_3d_net_traces(verts, edges)
            out_file = embedded_out_dir / f"{net_id_str}.json"
            with open(out_file, "w") as f:
                json.dump(traces, f)
            available_3d_ids.add(int(net_id_str))
        except Exception as e:
            print(f"Warning: Failed to generate 3D data for net {net_id_str}: {e}")

# 3. Generate parsed competitor JSON meshes in website/competitors/<net_id>/<stem>.json
site_competitors_dir = site_dir / "competitors"
competitor_json_count = 0

for net_id, comp_info in competitors_map.items():
    out_net_comp_dir = site_competitors_dir / str(net_id)
    out_net_comp_dir.mkdir(parents=True, exist_ok=True)

    json_files = []
    for txt_fname in comp_info["files"]:
        src_txt = data_dir / "competitors" / str(net_id) / txt_fname
        try:
            mesh_data = parse_competitor_file(src_txt)
            json_fname = f"{pathlib.Path(txt_fname).stem}.json"
            with open(out_net_comp_dir / json_fname, "w") as f:
                json.dump(mesh_data, f)
            json_files.append(json_fname)
            competitor_json_count += 1
        except Exception as e:
            print(f"Warning: Failed to parse competitor {src_txt}: {e}")

    comp_info["files"] = json_files

metadata = {
    "nets": nets,
    "competitors": competitors_map,
    "has_3d": sorted(list(available_3d_ids)),
}

output_json_path = site_dir / "site_metadata.json"
with open(output_json_path, "w") as f:
    json.dump(metadata, f)

print(
    f"Successfully processed {len(nets)} nets, "
    f"{len(available_3d_ids)} 3D embedded nets, "
    f"and {competitor_json_count} competitor meshes."
)

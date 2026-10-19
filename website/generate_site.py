import pathlib
import json
import sys

site_dir = pathlib.Path(__file__).resolve().parent
data_dir = site_dir.parent / "data"
scripts_dir = site_dir.parent / "scripts"

if str(scripts_dir) not in sys.path:
    sys.path.append(str(scripts_dir))

from visualize_spherical_geodesics import parse_metrization_file, build_3d_net_traces


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
                clean = fname.replace('_competitor.txt', '')
                parts = clean.split('_')
                assert int(parts[0]) == net_id
                if net_id == 827 and len(parts) != 5:
                    continue
                try:
                    [float(p) for p in parts[1:]]
                    valid_files.append(fname)
                except ValueError:
                    continue

            if valid_files:
                param_count = (len(valid_files[0].split('_')) - 2) // 2
                competitors_map[net_id] = {
                    "files": valid_files,
                    "param_count": param_count
                }

        # Status assignment precedence
        if net_id in special_red_ids:
            nets[net_id] = "red"
        elif has_metrization or (net_id in competitors_map):
            nets[net_id] = "green"
        else:
            nets[net_id] = "gray"

    return nets, competitors_map, embeddable_ids


nets, competitors_map, embeddable_ids = get_net_statuses()

# Generate 3D JSON trace data for embeddable nets
nets_3d_dir = site_dir / "net_3d_data"
nets_3d_dir.mkdir(parents=True, exist_ok=True)

available_3d_ids = set()
for net_id_str in embeddable_ids:
    metrized_txt = data_dir / "metrized_nets" / f"{net_id_str}_visualization.txt"
    if metrized_txt.exists():
        try:
            verts, edges = parse_metrization_file(metrized_txt)
            traces = build_3d_net_traces(verts, edges)
            out_file = nets_3d_dir / f"{net_id_str}.json"
            with open(out_file, "w") as f:
                json.dump(traces, f)
            available_3d_ids.add(int(net_id_str))
        except Exception as e:
            print(f"Warning: Failed to generate 3D data for net {net_id_str}: {e}")

metadata = {
    "nets": nets,
    "competitors": competitors_map,
    "has_3d": sorted(list(available_3d_ids))
}

output_json_path = site_dir / "site_metadata.json"
with open(output_json_path, "w") as f:
    json.dump(metadata, f)

print(f"Successfully processed {len(nets)} nets and generated {len(available_3d_ids)} 3D visual data files.")
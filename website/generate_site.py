import pathlib
import json

site_dir = pathlib.Path(__file__).resolve().parent
data_dir = site_dir.parent / "data"

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
                # Filter out the count header if present (e.g., '26')
                tokens = [t for t in content if t.isdigit()]
                # If the first token equals the count of remaining tokens, drop it
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
    special_red_ids = {"158", "681", "827", "920", "1029", "1100"}

    # 3. Evaluate all net IDs explicitly
    for net_int in range(total_nets):
        net_id = str(net_int)

        has_metrization = (net_id in embeddable_ids) or (data_dir / f"metrized_nets/{net_id}_visualization.txt").exists()

        net_comp_dir = competitors_dir / net_id
        if net_comp_dir.exists() and net_comp_dir.is_dir():
            all_files = [p.name for p in net_comp_dir.glob("*_competitor.txt")]

            # Filter valid standard competitor files and ignore legacy non-float files
            valid_files = []
            for fname in all_files:
                clean = fname.replace('_competitor.txt', '')
                parts = clean.split('_')
                try:
                    [float(p) for p in parts[1:]]
                    valid_files.append(fname)
                except ValueError:
                    continue

            if valid_files:
                # Hardcode parameter count or infer it
                if net_id == "827":
                    param_count = 2
                else:
                    sample_parts = valid_files[0].replace('_competitor.txt', '').split('_')[1:]
                    if len(sample_parts) >= 4:
                        param_count = 2
                    else:
                        param_count = 1 if len(sample_parts) == 1 else 0

                competitors_map[net_id] = {
                    "files": valid_files,
                    "param_count": param_count
                }

        # Status assignment precedence:
        if net_id in special_red_ids:
            nets[net_id] = "red"
        elif has_metrization or (net_id in competitors_map):
            nets[net_id] = "green"
        else:
            nets[net_id] = "gray"

    return nets, competitors_map

nets, competitors_map = get_net_statuses()

metadata = {
    "nets": nets,
    "competitors": competitors_map
}

output_json_path = site_dir / "site_metadata.json"
with open(output_json_path, "w") as f:
    json.dump(metadata, f)

print(f"Successfully processed {len(nets)} nets. Red count: {sum(1 for v in nets.values() if v == 'red')}")
import subprocess
import re

KMULT_PATH = "/home/smith/Compiled/muxmcmtool/kmult"
INPUT_FILE = "selections/2ADD/16on12bitsR100.txt"
OUTPUT_FILE = "kmult/16on12bitsR100.txt"
DOT_FILE = "module.dot"

def parse_dot_file():
    node_labels = {}
    edges = []

    with open(DOT_FILE, "r") as f:
        for line in f:
            line = line.strip()
            node_match = re.match(r'(node\d+)\s*\[.*label\s*=\s*"([^"]+)"', line)
            if node_match:
                node_id = node_match.group(1)
                label = node_match.group(2)
                node_labels[node_id] = label

            edge_match = re.match(r'(node\d+)\s*->\s*(node\d+)', line)
            if edge_match:
                src = edge_match.group(1)
                dst = edge_match.group(2)
                edges.append((src, dst))

    adder_count = sum(1 for label in node_labels.values() if label == "Add" or label == "Sub" or label == "AddSub")

    mux_incoming_sum = 0
    addsub_count = sum(1 for label in node_labels.values() if label == "AddSub")
    mux_incoming_sum += addsub_count
    
    mux_entries = {}
    for src, dst in edges:
        if dst in node_labels and node_labels[dst] == "Mux":
            mux_entries[dst] = mux_entries.get(dst, 0) + 1

    for entries in mux_entries.values():
        print(f"Mux with {entries} entries")
        mux_incoming_sum += entries - 1

    return adder_count, mux_incoming_sum

def main():
    total_adders = 0
    total_mux_entries = 0
    with open(INPUT_FILE, "r") as infile, open(OUTPUT_FILE, "w") as outfile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            nums = [s.strip() for s in line.split(";") if s.strip()]
            args = [KMULT_PATH] + ["-r100000 -c12 -f0"] + nums
            #args = [KMULT_PATH] + ["-c5 -f0"] + nums
            print(f"Running kmult with arguments: {args}")
            subprocess.run(args, check=True)
            
            try:
                adder_count, mux_incoming_sum = parse_dot_file()
                total_adders += adder_count
                total_mux_entries += mux_incoming_sum
                outfile.write(f"{adder_count} {mux_incoming_sum}\n")
                print(f"Processed line with numbers: {nums} -> Adders: {adder_count}, Mux entries sum: {mux_incoming_sum}")
            except Exception as e:
                print(f"Error processing DOT file for line '{line}': {e}")
        outfile.write(f"Average adders: {total_adders / 100}, Average mux entries: {total_mux_entries / 100}\n")

if __name__ == "__main__":
    main()

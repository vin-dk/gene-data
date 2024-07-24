import re

# file paths
mapping_file_path = r"C:\Users\13046\Desktop\gene_mapping.txt"
data_file_path = r"C:\Users\13046\Desktop\Work Paper\trials\steps\Trial 1_1\david_file_1.txt"
output_file_path = 'processed_data.txt'


def load_mapping(file_path):
    mapping = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            probe_id = parts[0]
            entrez_ids = [eid.strip() for eid in parts[1].split('///')]
            mapping[probe_id] = entrez_ids
    return mapping


def process_data(data_path, mapping):
    blocks = []
    with open(data_path, 'r') as f:
        block_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('Block'):
                if block_lines:
                    blocks.append('\n'.join(block_lines))
                    block_lines = []
                block_lines.append(line)
            elif line.startswith('Block Tqi:'):
                continue
            else:
                probe_ids = line.split(',')
                entrez_ids = []
                for probe_id in probe_ids:
                    probe_id = probe_id.strip()
                    if probe_id in mapping:
                        entrez_ids.extend(mapping[probe_id])
                # Remove duplicates and sort if needed
                entrez_ids = sorted(set(entrez_ids))
                block_lines.append(', '.join(entrez_ids))
        if block_lines:
            blocks.append('\n'.join(block_lines))
    return blocks


def write_output(output_path, blocks):
    with open(output_path, 'w') as f:
        for block in blocks:
            f.write(block + '\n\n')

mapping = load_mapping(mapping_file_path)
blocks = process_data(data_file_path, mapping)
write_output(output_file_path, blocks)

print(f"Processing complete. Output saved to {output_file_path}.")
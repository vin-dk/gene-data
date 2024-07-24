import csv

# File paths
input_file_path = r"C:\Users\13046\Desktop\GPL570-55999.txt"
output_file_path = r"C:\Users\13046\Desktop\gene_mapping.txt"


with open(input_file_path, 'r') as infile, open(output_file_path, 'w', newline='') as outfile:
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile)

   
    writer.writerow(['ID', 'Entrez'])
    
   
    for line in infile:
        if line.startswith('#'):
            continue
        infile.seek(0)
        break

    
    for row in reader:
        if len(row) >= 11:  
            probe_id = row[0]
            gene_symbol = row[11]

            
            writer.writerow([probe_id, gene_symbol])
        else:
            print(f"Skipped row with insufficient columns: {row}")
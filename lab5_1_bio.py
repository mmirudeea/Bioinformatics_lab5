def read_ncbi_fasta(file_path):
    sequence=[]
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('>'):
                sequence.append(line)
    return ''.join(sequence)

        

def brute_force_repeat_detector(sequence):
   
    seq_len = len(sequence)
    found_repeats = []

    for motif_len in range(3, 7):
        i = 0
        while i <= seq_len - motif_len:
            motif = sequence[i:i + motif_len]
            
            repeat_count = 1
            current_idx = i + motif_len
            
           
            while current_idx <= seq_len - motif_len:
                next_segment = sequence[current_idx : current_idx + motif_len]
                if next_segment == motif:
                    repeat_count += 1
                    current_idx += motif_len
                else:
                    break
            
            total_repetition_length = repeat_count * motif_len
            
           
            if repeat_count > 1 and 6 <= total_repetition_length <= 100:
                found_repeats.append({
                    'motif': motif,
                    'start_index': i,
                    'count': repeat_count,
                    'total_length': total_repetition_length
                })
                
                i = current_idx
            else:
                i += 1
                
    return found_repeats

if __name__ == "__main__":
    file_name = "gene.fna"
    dna_sequence = read_ncbi_fasta(file_name)
    
    if dna_sequence:
        print(f"Loaded DNA sequence of length {len(dna_sequence)} from {file_name}.")
        print("Running brute force engine to detect repetitions...\n")
        
        results = brute_force_repeat_detector(dna_sequence)
        
        if not results:
            print("No repetitions matching the criteria were found.")
        else:
            print(f"Found {len(results)} valid repetition blocks:\n")
            print(f"{'Motif':<10} | {'Start Idx':<10} | {'Repeats':<8} | {'Total Length':<12}")
            print("-" * 50)
            

            results_sorted = sorted(results, key=lambda x: x['start_index'])
            
            for r in results_sorted:
                print(f"{r['motif']:<10} | {r['start_index']:<10} | {r['count']:<8} | {r['total_length']:<12}")
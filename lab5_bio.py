import random 
import numpy as np
import matplotlib.pyplot as plt

def read_ncbi_fasta(file_path):
    sequence=[]
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('>'):
                sequence.append(line)
    return ''.join(sequence)

file_name="gene.fna"
dna_sequence = read_ncbi_fasta(file_name)
sequence_length=len(dna_sequence)

print(f"Read DNA sequence of {sequence_length}")

samples= []
num_samples=10

for _ in range (num_samples):
    max_frag_length=min(3000, sequence_length)
    if sequence_length < 100:
        frag_length = sequence_length
    else:
        frag_length = random.randint(100, max_frag_length)
        
    start_idx = random.randint(0, sequence_length-frag_length)
    fragment = dna_sequence[start_idx : start_idx +frag_length]
    samples.append(fragment)
    
fragment_lengths = [len(frag) for frag in samples]
print(f"Lengths of the 10 random fragments: {fragment_lengths}")


fig, ax = plt.subplots(figsize=(5,8))

max_expected_len = 3000
min_expected_len = 100

migration_distances = [np.log10(max_expected_len) - np.log10(length) for length in fragment_lengths if length > 0]

lane_center = 0.5
band_width = 0.6

for dist in migration_distances:
    ax.plot(
        [lane_center - band_width/2, lane_center + band_width/2],
        [dist, dist],
        color = '#1d2024',
        linewidth=4,
        alpha=0.8
        )

ax.set_xlim(0,1)

max_travel = np.log10(max_expected_len) - np.log10(min_expected_len)
ax.set_ylim(-0.1, max_travel + 0.2)
ax.invert_yaxis() 
ax.set_xticks([])

ladder_lengths = [3000, 2000, 1500, 1000, 500, 250, 100]
ladder_distances = [np.log10(max_expected_len) - np.log10(l) for l in ladder_lengths]
ax.set_yticks(ladder_distances)
ax.set_yticklabels([f"{l} bp" for l in ladder_lengths])

ax.set_ylabel("Fragment Size (Base Pairs)")
ax.set_title("Simulated Gel Electrophoresis\n(Migration based on molecular weight)")

ax.plot(
    [lane_center - band_width/2, lane_center + band_width/2], 
    [0, 0], 
    color='gray', 
    linewidth=12, 
    solid_capstyle='butt'
)
ax.text(lane_center, -0.06, 'Injection Well', ha='center', va='bottom', fontweight='bold')


ax.set_facecolor('#e8f4f8') 
plt.tight_layout()
plt.show()
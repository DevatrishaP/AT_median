#Make combined fasta file for genes whose median AT content is to be plotted

# Install necessary packages
# !pip install biopython
# !pip install scipy

from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Parameters
bin_size = 500  # Number of bins to split each gene sequence into equal parts
output_csv = 'combined_at_content.csv'
output_plot = 'combined_at_content_plot.png'
fasta_files = {
    'File1': 'path_file1', 
    'File2': 'path_file2',
    'File3': 'path_file3'
}
#num_interpolated_points = 100  # Number of points for smoothing ##use only if you want the lines in the plot to be smoothened out

# Function to calculate AT content for each fixed bin in a gene sequence
def calculate_at_content_fixed_bins(fasta_file, label):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    data = []

    # Loop through each sequence in the FASTA file
    for record in sequences:
        gene_id = record.id
        seq = str(record.seq)
        seq_length = len(seq)

        # Divide sequence into equal bins
        bin_indices = np.linspace(0, seq_length, num=bin_size + 1, dtype=int)

        # Calculate AT content for each bin
        for i in range(len(bin_indices) - 1):
            start, end = bin_indices[i], bin_indices[i + 1]
            bin_seq = seq[start:end]
            
            # Calculate AT content for the bin
            at_count = sum(1 for base in bin_seq if base == "A" or base == "T")
            at_content = at_count / len(bin_seq) if len(bin_seq) > 0 else 0  # Avoid division by zero

            # Append data
            data.append([label, gene_id, i + 1, at_content])  # Record file label, gene_id, bin_index, and AT content

    # Create DataFrame and return
    df = pd.DataFrame(data, columns=['File', 'Gene', 'Bin', 'AT_Content'])
    return df

# Calculate AT content for each file and combine results
df_file1 = calculate_at_content_fixed_bins(fasta_files['File1'], 'File1')
df_file2 = calculate_at_content_fixed_bins(fasta_files['File2'], 'File2')
df_file3 = calculate_at_content_fixed_bins(fasta_files['File3'], 'File3')
df_combined = pd.concat([df_file1, df_file2, df_file3])

# Save combined data to CSV
df_combined.to_csv(output_csv, index=False)
print(f"Data saved to {output_csv}")

# Plotting only the smoothed median AT content for each file
plt.figure(figsize=(10, 6))

# Function to smooth a line using interpolation
def smooth_line(data, color, label):
    x = data.index
    y = data.values
    interpolator = interp1d(x, y, kind='cubic')
    x_smooth = np.linspace(x.min(), x.max(), num_interpolated_points)
    y_smooth = interpolator(x_smooth)
    plt.plot(x_smooth, y_smooth, color=color, linewidth=1.5, label=label)

# Calculate and plot smoothed median AT content for each file
median_file1 = df_file1.groupby('Bin')['AT_Content'].median()
#smooth_line(median_file1, 'red', "Median AT Content (file1)") ##use only when doing interpolation
plt.plot(median_file1.index, median_file1.values, color='red', linewidth=1.5, label="Median AT Content (File1)")

median_file2 = df_file2.groupby('Bin')['AT_Content'].median()
#smooth_line(median_file2, 'blue', "Median AT Content (file2)") ##use only when doing interpolation
plt.plot(median_file2.index, median_file2.values, color='blue', linewidth=1.5, label="Median AT Content (File2)")

median_file3 = df_file3.groupby('Bin')['AT_Content'].median()
#smooth_line(median_file3, 'green', "Median AT Content (file3)") ##use only when doing interpolation
plt.plot(median_file3.index, median_file3.values, color='green', linewidth=1.5, label="Median AT Content (File3)")

# Customize plot
plt.xlabel('Bin')
plt.ylabel('AT Content')
plt.title('Smoothed Median AT Content Across Genes')
plt.legend()
plt.savefig(output_plot)
plt.show()

print(f"Plot saved as {output_plot}")

import matplotlib.pyplot as plt
import pandas as pd

# Load data from TSV file
data = pd.read_csv('1nql_biopython.tsv', delimiter='\t',
                   header=None, names=['AA', 'Phi', 'Psi', 'Type'])

# Create Ramachandran plot using Matplotlib
plt.scatter(data['Phi'], data['Psi'], s=10, c='black')
plt.xlabel('Phi')
plt.ylabel('Psi')
plt.title('Ramachandran Plot for Protein 1NQL')
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.xticks([-180, -120, -60, 0, 60, 120, 180])
plt.yticks([-180, -120, -60, 0, 60, 120, 180])
plt.grid(True)
plt.show()

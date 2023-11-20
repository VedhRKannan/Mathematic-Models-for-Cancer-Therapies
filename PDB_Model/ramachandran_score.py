# Import the csv module to read tab-separated values (TSV) files.
import csv

# Define a function that reads phi and psi angles from a TSV file.


def read_phi_psi_angles(file_name):
    angles = []
    # Open the TSV file in read mode.
    with open(file_name, 'r') as file:
        # Create a csv reader object and specify the delimiter as tab.
        reader = csv.reader(file, delimiter='\t')
        # Loop through each row in the file.
        for row in reader:
            # Extract the phi and psi angles from the second and third columns, respectively.
            phi, psi = float(row[1]), float(row[2])
            # Add the phi and psi angles as a tuple to the angles list.
            angles.append((phi, psi))
    # Return the angles list.
    return angles

# Define a function that counts the number of phi/psi angle pairs in each region of the Ramachandran plot.


def count_regions(angles):
    # Initialize the counts for the allowed, favored, and disallowed regions.
    allowed = 0
    favored = 0
    disallowed = 0

    # Loop through each phi/psi angle pair in the angles list.
    for phi, psi in angles:
        # Check if the angle pair falls within the favored region.
        if (-135 <= phi <= -45) and (-90 <= psi <= 175):
            favored += 1
        # Check if the angle pair falls within the allowed region.
        elif (-180 <= phi <= 180) and (-180 <= psi <= 180):
            allowed += 1
        # Otherwise, the angle pair falls within the disallowed region.
        else:
            disallowed += 1

    # Return the counts for each region.
    return allowed, favored, disallowed

# Define a function that calculates the Ramachandran plot quality (RPQ) for a given TSV file.


def calculate_rpq(file_name):
    # Read the phi/psi angles from the TSV file.
    angles = read_phi_psi_angles(file_name)
    # Count the number of angle pairs in each region of the Ramachandran plot.
    allowed, favored, disallowed = count_regions(angles)
    # Calculate the total number of residues.
    total_residues = len(angles)
    # Calculate the percentage of residues in the allowed, favored, and disallowed
    allowed_percent = (allowed / total_residues) * 100
    favored_percent = (favored / total_residues) * 100
    disallowed_percent = (disallowed / total_residues) * 100

    rpq = favored_percent + 0.5 * allowed_percent

    return rpq, allowed_percent, favored_percent, disallowed_percent


if __name__ == "__main__":
    file_name = "Alphafold/prediction_biopython.tsv"
    rpq, allowed, favored, disallowed = calculate_rpq(file_name)

    print(f"RPQ: {rpq:.2f}%")
    print(f"Allowed region: {allowed:.2f}%")
    print(f"Favored region: {favored:.2f}%")
    print(f"Disallowed region: {disallowed:.2f}%")

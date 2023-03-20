import csv


def read_phi_psi_angles(file_name):
    angles = []
    with open(file_name, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            phi, psi = float(row[1]), float(row[2])
            angles.append((phi, psi))
    return angles


def count_regions(angles):
    allowed = 0
    favored = 0
    disallowed = 0

    for phi, psi in angles:
        if (-135 <= phi <= -45) and (-90 <= psi <= 175):
            favored += 1
        elif (-180 <= phi <= 180) and (-180 <= psi <= 180):
            allowed += 1
        else:
            disallowed += 1

    return allowed, favored, disallowed


def calculate_rpq(file_name):
    angles = read_phi_psi_angles(file_name)
    allowed, favored, disallowed = count_regions(angles)
    total_residues = len(angles)

    allowed_percent = (allowed / total_residues) * 100
    favored_percent = (favored / total_residues) * 100
    disallowed_percent = (disallowed / total_residues) * 100

    rpq = favored_percent + 0.5 * allowed_percent

    return rpq, allowed_percent, favored_percent, disallowed_percent


if __name__ == "__main__":
    file_name = "new.tsv"
    rpq, allowed, favored, disallowed = calculate_rpq(file_name)

    print(f"RPQ: {rpq:.2f}%")
    print(f"Allowed region: {allowed:.2f}%")
    print(f"Favored region: {favored:.2f}%")
    print(f"Disallowed region: {disallowed:.2f}%")

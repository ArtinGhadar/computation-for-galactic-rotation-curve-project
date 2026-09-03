import csv

with open("sparc_data_dmf.csv", "r", newline="") as file:
    reader = csv.reader(
        line for line in file
        if line.strip() and not line.lstrip().startswith("#")
    )

    database = list(reader)

galaxy_inp = 'D512-2'
variable_inp = 'Vobs'
variable_map = {
    'Radius': 0,
    'Vobs': 1,
    'errV': 2,
    'Vgas': 3,
    'Vdisk': 4,
    'Vbul': 5,
    'SBdisk': 6,
    'SBbul': 7
}
counter = -1
data = []
for row in database:
    counter += 1
    if row[0] == galaxy_inp:
        data.append(float(row[variable_map[variable_inp]]))
print(data)
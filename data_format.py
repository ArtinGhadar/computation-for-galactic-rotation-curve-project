import csv

variable_map = {
        'Galaxy': 0,
        'Radius': 1,
        'Vobs': 2,
        'errV': 3,
        'Vgas': 4,
        'Vdisk': 5,
        'Vbul': 6,
        'SBdisk': 7,
        'SBbul': 8
    }

with open("sparc_data_dmf.csv", "r", newline="") as file:
    reader = csv.reader(
        line for line in file
        if line.strip() and not line.lstrip().startswith("#")
    )

    database = list(reader)

def data_format(galaxy_inp, variable_inp):
    
    data = []
    for row in database: 
        if row[0] == galaxy_inp:
            data.append(float(row[variable_map[variable_inp]]))

    return data
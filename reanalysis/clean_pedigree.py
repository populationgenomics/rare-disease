"""
read a ped file, write without any family structure
removes associations between parents and children
re-assigns a personal family ID to each participant
"""

PED = 'input/acute_care_updated.ped'

output_lines = []

with open(PED, 'r', encoding='utf-8') as read_handle:

    for line in read_handle:

        if line.startswith('#'):
            output_lines.append(line)
            continue

        # split the string into a list
        llist = line.split('\t')

        # treat the individual as a 'family'
        llist[0] = llist[1]

        # delete parents
        llist[2] = ''
        llist[3] = ''
        output_lines.append('\t'.join(llist))

with open('input/no_family_structure.ped', 'w', encoding='utf-8') as handle:
    handle.writelines(output_lines)

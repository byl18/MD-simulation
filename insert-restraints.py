import sys
import shutil

restraints_directives = {}

for i in range(1, len(sys.argv)):
    if i % 2 == 1:
        restraints_directives[int(sys.argv[i])] = sys.argv[i + 1]

shutil.copy("gleap.top", "gleap.before_insert_restraints.top")

new_top = ""
molecule_type_count = 0
for line in open("gleap.top").readlines():
    if line.strip() == "[ moleculetype ]":
        molecule_type_count += 1

    if molecule_type_count in restraints_directives:
        new_top += (
            "\n\n" + open(restraints_directives[molecule_type_count]).read() + "\n\n"
        )
        del restraints_directives[molecule_type_count]

    new_top += line
open("gleap.top", "w").write(new_top)

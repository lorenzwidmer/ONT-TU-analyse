
import re
import csv

bsaI_regex = re.compile('GGTCTC[ACGT]([ACGT]{4})([ACGT]+)([ACGT]{4})[ACGT]GAGACC', flags=re.IGNORECASE)
lvl0 = ["CTAT", "ACTT", "CATA", "GGAA", "AGTG", "ACCG", "GGCT", "CGAC", "TGTT"]
lvl0_dict = {overhang:i for i,overhang in enumerate(lvl0)}
partnames = ['promo', 'rbs', 'signal', 'ntag', 'cds', 'ctag', 'stop', 'term']
output = []

with open('mocolo_parts.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='"')
    for row in spamreader:
        fullname, author, seq, lvl = row

        if lvl != 'lvl0':
            continue

        matsch = bsaI_regex.search(seq.upper())
        if matsch is not None:
            start = lvl0_dict.get(matsch.group(1), 0)
            end = lvl0_dict.get(matsch.group(3), 0)

            name = fullname.strip('✅☑️ ').replace(' ', '-')
            parts = '-'.join(partnames[start:end])
            print(parts)

            output.append(f">{name}|{parts}\n")
            output.append(f"{matsch.group(2)}\n")
        else:
            print(f'{fullname} contains no bsai')

with open('parts.fasta', 'w') as f:
    f.writelines(output)

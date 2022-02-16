from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

s1 = SeqRecord(id="AID", seq=Seq("ACGT"), description="")

from random import choice, choices, randint
from input_parser import InputParser
import re

def msa_neighbor_add_remove(df, changes=1):
    seq_count = len(df.index)
    if changes > seq_count:
        changes = seq_count

    random_seq_indexes = choices(list(df.index), k=changes)

    with open("/tmp/msa_neighbor_add_remove.fasta", "w") as tmp_file:
        dft = df.transpose()
        for index, row in dft.items():
            values = dft[index]
            sequence = ''.join([x for x in values])

            if index in random_seq_indexes:
                if randint(0, 1) == 0:
                    pos = randint(0, len(sequence) - 1)
                    new_sequence = "".join((sequence[:pos], "-", sequence[pos:]))
                else:
                    gap_positions = [i.start() for i in re.finditer("-", sequence)]
                    if len(gap_positions) > 0:
                        pos = choice(gap_positions)
                        new_sequence = sequence[0:pos:] + sequence[pos + 1::]
                    else:
                        new_sequence = sequence

                sequence = new_sequence

            tmp_file.write(">{0}\n{1}\n".format(index, sequence))

    sequences_dictionary = InputParser.read_fasta_to_dict(["/tmp/msa_neighbor_add_remove.fasta"])
    neighbor = InputParser.build_dataframe(sequences_dictionary)

    return neighbor

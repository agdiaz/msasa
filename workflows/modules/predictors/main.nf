process predictMuscle {
    tag "$iteration"
    debug true

    input:
    tuple val(iteration), path(sequences)
    
    output:
    path "*_muscle.fasta", emit: alignment
    path "*_muscle.time",  emit: timing

    script:
    """
    /usr/bin/time -p -o ${sequences.simpleName}_muscle.time /software/muscle -align ${sequences} -amino -output output.fasta
    fold -w 60 output.fasta > ${sequences.simpleName}_muscle.fasta
    """ 
}

process sortPredictedAlignment {
    tag "${referenceAlignment.simpleName} -> ${disorderedAlignment.simpleName}"
    
    input:
    path referenceAlignment
    path disorderedAlignment

    output:
    path "*_sorted.fasta", emit: alignment

    script:
    """
    #!/usr/local/bin/python

    class FastaSequence:
        def __init__(self, id, sequence):
            self.id = id
            self.sequence = sequence

        def __repr__(self):
            return '>' + self.id + '\\n' + self.sequence

    class FastaFile:
        def __init__(self, filename):
            self.filename = filename
            self.sequences = []

        def read(self):
            print('Reading', self.filename)
            with open(self.filename, 'r') as file:
                seq = ''
                name = ''
                for line in file:
                    line = line.strip()
                    if line.startswith('>'):
                        if name:
                            self.sequences.append(FastaSequence(name, seq))
                        name = line[1:]
                        seq = ''
                    else:
                        seq += line
                self.sequences.append(FastaSequence(name, seq))

        def write(self, sequences):
            print('Writing', self.filename)

            with open(self.filename, 'w') as file:
                for sequence in sequences:
                    file.write(str(sequence) + '\\n')

    # Read the sequences from file1.fasta
    file1 = FastaFile('${disorderedAlignment}')
    file1.read()

    # Read the sequence IDs from file2.fasta
    sequence_ids = []
    with open('${referenceAlignment}', 'r') as file:
        for line in file:
            if line.startswith('>'):
                sequence_ids.append(line[1:].strip())

    # Sort the sequences from file1.fasta based on the order of the sequence IDs from file2.fasta
    print('sequence_ids', sequence_ids)
    sorted_sequences = [sequence for id in sequence_ids for sequence in file1.sequences if sequence.id == id]
    print('sorted_sequences', sorted_sequences)

    # Write the sorted sequences to a new file
    file2 = FastaFile('${disorderedAlignment.simpleName}_sorted.fasta')
    file2.write(sorted_sequences)
    """
}

process predictClustalOmega {
    tag "$iteration"

    input:
    tuple val(iteration), path(sequences)
    
    output:
    path "*_clustalo.fasta", emit: alignment
    path "*_clustalo.time",  emit: timing

    script:
    """
    #!/bin/bash
    /usr/bin/time -p -o ${sequences.simpleName}_clustalo.time /software/clustalo --force -i ${sequences} -o ${sequences.simpleName}_clustalo.fasta
    """ 
}

process predictMafft {
    tag "$iteration"

    input:
    tuple val(iteration), path(sequences)
    
    output:
    path "*_mafft.fasta", emit: alignment
    path "*_mafft.time",  emit: timing

    script:
    """
    #!/bin/bash
    
    /usr/bin/time -p -o ${sequences.simpleName}_mafft.time /usr/local/bin/mafft --auto ${sequences} > ${sequences.simpleName}_mafft.fasta
    """ 
}

process predictTCoffee {
    tag "$iteration"

    input:
    tuple val(iteration), path(sequences)
    
    output:
    path "*_tcoffee.fasta", emit: alignment
    path "*_tcoffee.time",  emit: timing

    script:
    """
    t_coffee -seq ${sequences}
    /usr/bin/time -p --output=${sequences.simpleName}_tcoffee.time t_coffee -other_pg seq_reformat -in ${sequences.simpleName}.aln -output fasta_aln > ${sequences.simpleName}_tcoffee.fasta
    """
}

process predictHmmer {
    tag "$iteration - $referenceFasta"

    input:
    tuple val(iteration), path(sequences)
    path referenceFasta
    
    output:
    path "*_hmmer.fasta", emit: alignment
    path "*_hmmer.time",  emit: timing

    script:
    """
    /usr/bin/time -p -o ${sequences.simpleName}_hmmer.time cp ${referenceFasta} ${sequences.simpleName}_hmmer.fasta

    # hmmbuild --amino reference.hmm ${referenceFasta}
    # /usr/bin/time -p -o ${sequences.simpleName}_hmmer.time hmmalign --amino reference.hmm ${sequences} > ${sequences.simpleName}_hmmer.sto
    # esl-reformat fasta ${sequences.simpleName}_hmmer.sto > ${sequences.simpleName}_hmmer.fasta
    """
}

process predictKAlign {
    tag "$iteration"

    input:
    tuple val(iteration), path(sequences)
    
    output:
    path "*_kalign.fasta", emit: alignment
    path "*_kalign.time",  emit: timing

    script:
    """
    #!/bin/bash
    /usr/bin/time -p -o ${sequences.simpleName}_kalign.time /software/kalign -input ${sequences} -output ${sequences.simpleName}_kalign.fasta -format fasta -quiet
    """
}

process predictMsasaSingleMs {
    errorStrategy 'ignore'

    input:
    tuple val(iteration), path(sequences)

    output:
    path "*_msasa_singlems.fasta", emit: alignment
    path "*_msasa_singlems.time",  emit: timing

    script:
    """
    /usr/local/bin/python /usr/src/app/src/msa.py --input ${sequences} \
        --output output.fasta \
        --comparer single_ms \
        --n-iterations 100 \
        --temperature 10 \
        --execution-id 1 > ${sequences.simpleName}_msasa_singlems.time
    
    fold -w 60 output.fasta > ${sequences.simpleName}_msasa_singlems.fasta
    """
}

process predictMsasaBlosum {
    errorStrategy 'ignore'
    
    input:
    tuple val(iteration), path(sequences)

    output:
    path "*_msasa_singleblosum.fasta", emit: alignment
    path "*_msasa_singleblosum.time",  emit: timing

    script:
    """
    /usr/local/bin/python /usr/src/app/src/msa.py --input ${sequences} \
        --output output.fasta \
        --comparer single_blosum \
        --n-iterations 100 \
        --temperature 20 \
        --execution-id 1 > ${sequences.simpleName}_msasa_singleblosum.time
    
    fold -w 60 output.fasta > ${sequences.simpleName}_msasa_singleblosum.fasta
    """
}

process predictMsasaSingleMatching {
    errorStrategy 'ignore'

    input:
    tuple val(iteration), path(sequences)

    output:
    path "*_msasa_singlematching.fasta", emit: alignment
    path "*_msasa_singlematching.time",  emit: timing

    script:
    """
    /usr/local/bin/python /usr/src/app/src/msa.py --input ${sequences} \
        --output output.fasta \
        --comparer single_matching \
        --n-iterations 100 \
        --temperature 20 \
        --execution-id 1 > ${sequences.simpleName}_msasa_singlematching.time
    
    fold -w 60 output.fasta > ${sequences.simpleName}_msasa_singlematching.fasta
    """
}

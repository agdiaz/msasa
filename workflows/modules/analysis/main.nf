params.referenceXml = ""
reference = file(params.referenceXml)

process computeBaliScore {
    tag "${alignment}"
    errorStrategy 'ignore'

    publishDir "$baseDir/tesis_results/${reference.simpleName}", mode: 'copy'

    input:
        path alignment

    output:
        path "*.baliscore.txt", emit: score

    shell:
    '''
    /software/bali-score/target/release/bali-score -t !{alignment} -r !{reference} -o !{alignment.baseName}.baliscore.txt
    '''
}

process computeMumsaOverlapScore {
    tag "${referenceAlignment}"
    publishDir "$baseDir/tesis_results/${reference.simpleName}", mode: 'copy'

    input:
        path referenceAlignment
        path muscleAlignment
        path clustalAlignment
        path mafftAlignment
        path tCoffeeAlignment
        path hmmerAlignment
        path kalignAlignment

    output:
        path "${referenceAlignment.simpleName}_mumsa.txt"

    script:
    """
    /software/mumsa-1.0/mumsa -r -q $referenceAlignment \
        $muscleAlignment \
        $clustalAlignment \
        $mafftAlignment \
        $tCoffeeAlignment \
        $kalignAlignment \
        > ${referenceAlignment.simpleName}_mumsa.txt

    cat ${referenceAlignment.simpleName}_mumsa.txt
    """
}

process exportAlignments {
    tag "${reference.simpleName}"

    publishDir "$baseDir/tesis_results/${reference.simpleName}", mode: 'copy'

    input:
    path muscleAlignment
    path clustalAlignment
    path mafftAlignment
    path tCoffeeAlignment
    path hmmerAlignment
    path kalignAlignment

    output:
    path "*.fasta"

    script:
    """
    cp $muscleAlignment ${reference.simpleName}_muscle.fasta
    cp $clustalAlignment ${reference.simpleName}_clustal.fasta
    cp $mafftAlignment ${reference.simpleName}_mafft.fasta
    cp $tCoffeeAlignment ${reference.simpleName}_tCoffee.fasta
    cp $hmmerAlignment ${reference.simpleName}_hmmer.fasta
    cp $kalignAlignment ${reference.simpleName}_kalign.fasta
    """
}
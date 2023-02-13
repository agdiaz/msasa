params.referenceXml = ""
reference = file(params.referenceXml)

process computeCoreIndex {
    tag "${alignment.simpleName}"
    publishDir "$baseDir/tesis_results/${reference.simpleName}", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        val predictorName
        path alignment

    output:
        path "${reference.simpleName}_coreindex_${predictorName}.txt"

    script:
    """
    t_coffee -infile=$alignment -output=score_ascii -score -outfile=${reference.simpleName}_coreindex_${predictorName}.txt
    """
}

process computeBaliScore {
    tag "${alignment.simpleName}"
    publishDir "$baseDir/tesis_results/${reference.simpleName}", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        val predictorName
        path alignment

    output:
        path "${reference.simpleName}_baliscore_${predictorName}.txt", emit: score

    shell:
    """
    /software/bali-score/target/release/bali-score -t ${alignment} \
        -r ${reference} \
        -o ${reference.simpleName}_baliscore_${predictorName}.txt
    """
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
    path msasaSingleMsAlignment
    path msasaBlosumAlignment
    path msasaSingleMatchingAlignment

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
        $msasaSingleMsAlignment \
        $msasaBlosumAlignment \
        $msasaSingleMatchingAlignment \
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
    path msasaSingleMsAlignment
    path msasaBlosumAlignment
    path msasaSingleMatchingAlignment

    output:
    path "*.fasta"

    script:
    """
    ln -s $muscleAlignment ${reference.simpleName}_alignment_muscle.fasta
    ln -s $clustalAlignment ${reference.simpleName}_alignment_clustal.fasta
    ln -s $mafftAlignment ${reference.simpleName}_alignment_mafft.fasta
    ln -s $tCoffeeAlignment ${reference.simpleName}_alignment_tCoffee.fasta
    ln -s $hmmerAlignment ${reference.simpleName}_alignment_hmmer.fasta
    ln -s $kalignAlignment ${reference.simpleName}_alignment_kalign.fasta
    ln -s $msasaSingleMsAlignment ${reference.simpleName}_alignment_msasa_singlems.fasta
    ln -s $msasaBlosumAlignment ${reference.simpleName}_alignment_msasa_blosum.fasta
    ln -s $msasaSingleMatchingAlignment ${reference.simpleName}_alignment_msasa_singlematching.fasta
    """
}
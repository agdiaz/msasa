params.data = '/home/adrian/workspace/msasa/test/small_short.fa'

process runClustalOmega {
    publishDir "$baseDir/results", mode: 'copy'

    input:
        path inputFile

    output:
        path "*.clustalo.afa"

    shell:
    """
    /software/clustalo --force -i $inputFile -o result.clustalo.afa
    """
}

process runKalign {
    publishDir "$baseDir/results", mode: 'copy'

    input:
        path inputFile

    output:
        path "*.kalign.afa"

    shell:
    """
    /software/kalign -i $inputFile -o result.kalign.afa -format fasta
    """
}

process runMAFFT {
    publishDir "$baseDir/results", mode: 'copy'

    input:
        path inputFile

    output:
        path "*.mafft.afa"

    shell:
    """
    /usr/local/bin/mafft --auto $inputFile > result.mafft.afa
    """
}

process runMuscle {
    publishDir "$baseDir/results", mode: 'copy'

    input:
        path inputFile

    output:
        path "*.muscle.afa"

    shell:
    """
    /software/muscle -align $inputFile -output longresult.afa
    fold -w 60 -s longresult.afa > result.muscle.afa
    """
}

process runTCoffee {
    publishDir "$baseDir/results", mode: 'copy'

    input:
        path inputFile

    output:
        path "*.t_coffee.afa"

    shell:
    """
    t_coffee -seq $inputFile
    t_coffee -other_pg seq_reformat -in ${inputFile.baseName}.aln -output fasta_aln > ${inputFile.baseName}.t_coffee.afa
    """
}

process runMSASA {
    publishDir "$baseDir/results", mode: 'copy'

    input:
        path inputFile

    output:
        path "*.msasa-blosum.afa"
        path "*.msasa-globalms.afa"
        path "*.msasa-matching.afa"

    shell:
    """
    python3.10 /software/msasa/src/msa.py --input $inputFile \
        --output longresult-blosum.afa \
        --comparer blosum \
        --n-iterations 500 \
        --optimization min

    python3.10 /software/msasa/src/msa.py --input $inputFile \
        --output longresult-globalms.afa \
        --comparer global_ms_min \
        --n-iterations 50 \
        --optimization min

    python3.10 /software/msasa/src/msa.py --input $inputFile \
        --output longresult-matching.afa \
        --comparer matching \
        --n-iterations 7000 \
        --n-iterations 500 \
        --temperature 1 \
        --optimization min

    fold -w 60 -s longresult-blosum.afa > result.msasa-blosum.afa
    fold -w 60 -s longresult-globalms.afa > result.msasa-globalms.afa
    fold -w 60 -s longresult-matching.afa > result.msasa-matching.afa
    """
}


process compareWithBaseline {
    debug true
    publishDir "$baseDir/results", mode: 'copy'

    input:
        path clustalo
        path kalign
        path mafft
        path muscle
        path t_coffee
        path msasaBlosum
        path msasaGlobalMs
        path msasaMatching

    shell:
    """
    echo clustalo:
    cat $clustalo

    echo kalign:
    cat $kalign

    echo mafft:
    cat $mafft

    echo muscle:
    cat $muscle

    echo t_coffee:
    cat $t_coffee

    echo msasa-blosum:
    cat $msasaBlosum

    echo msasaGlobalMs:
    cat $msasaGlobalMs

    echo msasaMatching:
    cat $msasaMatching

    t_coffee -infile=$clustalo -output=html -score
    ls -lah .
    """
}

workflow {
    runClustalOmega(params.data)
    runKalign(params.data)
    runMAFFT(params.data)
    runMuscle(params.data)
    runTCoffee(params.data)
    runMSASA(params.data)

    compareWithBaseline(
        runClustalOmega.out,
        runKalign.out,
        runMAFFT.out,
        runMuscle.out,
        runTCoffee.out,
        runMSASA.out
    )
}

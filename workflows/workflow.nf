#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.data = '/home/adrian/workspace/msasa/test/small_short.fa'
params.outputDir = 'small_short'

process runClustalOmega {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

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
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

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
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

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
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

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
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

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
    debug true
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path inputFile

    output:
        path "*.msasa-blosum.afa"
        path "*.msasa-globalms.afa"
        path "*.msasa-matching.afa"

    shell:
    '''
    LOG_BLOSUM_FILENAME=!{inputFile.baseName}_blosum.log
    LOG_BLOSUM_FILENAME_SORTED=$LOG_BLOSUM_FILENAME.sorted

    LOG_GLOBALMS_FILENAME=!{inputFile.baseName}_globalms.log
    LOG_GLOBALMS_FILENAME_SORTED=$LOG_GLOBALMS_FILENAME.sorted

    LOG_MATCHING_FILENAME=!{inputFile.baseName}_matching.log
    LOG_MATCHING_FILENAME_SORTED=$LOG_MATCHING_FILENAME.sorted


    for ((i=1;i<=3;i++));
    do
        OUTPUT_MSA_BLOSUM_FILENAME=$i.msasa_blosum.afa
        OUTPUT_MSA_GLOBALMS_FILENAME=$i.msasa_globalms.afa
        OUTPUT_MSA_MATCHING_FILENAME=$i.msasa_matching.afa

        echo "MSASA Execution # $i started"
        python3.10 /software/msasa/src/msa.py --input !{inputFile} \
        --output $OUTPUT_MSA_BLOSUM_FILENAME \
        --comparer blosum \
        --n-iterations 500 \
        --optimization min >> $LOG_BLOSUM_FILENAME

        python3.10 /software/msasa/src/msa.py --input !{inputFile} \
        --output $OUTPUT_MSA_GLOBALMS_FILENAME \
        --comparer global_ms_min \
        --n-iterations 50 \
        --optimization min >> $LOG_GLOBALMS_FILENAME

        python3.10 /software/msasa/src/msa.py --input !{inputFile} \
        --output $OUTPUT_MSA_MATCHING_FILENAME \
        --comparer matching \
        --n-iterations 7000 \
        --n-iterations 500 \
        --temperature 1 \
        --optimization min >> $LOG_MATCHING_FILENAME
    done

    cat $LOG_BLOSUM_FILENAME | sort -r --field-separator=";" --key=3 > $LOG_BLOSUM_FILENAME_SORTED
    cat $LOG_GLOBALMS_FILENAME | sort -r --field-separator=";" --key=3 > $LOG_GLOBALMS_FILENAME_SORTED
    cat $LOG_MATCHING_FILENAME | sort -r --field-separator=";" --key=3 > $LOG_MATCHING_FILENAME_SORTED

    BEST_BLOSUM=$(cat $LOG_BLOSUM_FILENAME_SORTED | head -n 1 | cut -d ';' -f1 )
    BEST_GLOBALMS=$(cat $LOG_GLOBALMS_FILENAME_SORTED | head -n 1 | cut -d ';' -f1 )
    BEST_MATCHING=$(cat $LOG_MATCHING_FILENAME_SORTED | head -n 1 | cut -d ';' -f1 )

    echo BEST_BLOSUM: $BEST_BLOSUM
    cp $BEST_BLOSUM longresult-blosum.afa

    echo BEST_GLOBALMS: $BEST_GLOBALMS
    cp $BEST_GLOBALMS longresult-globalms.afa

    echo BEST_MATCHING: $BEST_MATCHING
    cp $BEST_MATCHING longresult-matching.afa

    fold -w 60 -s longresult-blosum.afa > result.msasa-blosum.afa
    fold -w 60 -s longresult-globalms.afa > result.msasa-globalms.afa
    fold -w 60 -s longresult-matching.afa > result.msasa-matching.afa
    '''
}


process computeCoreIndex {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path clustalo
        path kalign
        path mafft
        path muscle
        path t_coffee
        path msasaBlosum
        path msasaGlobalMs
        path msasaMatching

    output:
        path "*.core-index.html"

    shell:
    """
    t_coffee -infile=$clustalo -output=html -score -outfile=clustalo.core-index.html
    t_coffee -infile=$kalign -output=html -score -outfile=kalign.core-index.html
    t_coffee -infile=$mafft -output=html -score -outfile=mafft.core-index.html
    t_coffee -infile=$muscle -output=html -score -outfile=muscle.core-index.html
    t_coffee -infile=$t_coffee -output=html -score -outfile=t_coffee.core-index.html
    t_coffee -infile=$msasaBlosum -output=html -score -outfile=msasaBlosum.core-index.html
    t_coffee -infile=$msasaGlobalMs -output=html -score -outfile=msasaGlobalMs.core-index.html
    t_coffee -infile=$msasaMatching -output=html -score -outfile=msasaMatching.core-index.html
    """
}

process computeTransitiveConsistencyScore {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path clustalo
        path kalign
        path mafft
        path muscle
        path t_coffee
        path msasaBlosum
        path msasaGlobalMs
        path msasaMatching

    output:
        path "*.tcs.html"

    shell:
    """
    t_coffee -infile $clustalo -evaluate -output=score_html -outfile=clustalo.tcs.html
    t_coffee -infile $kalign -evaluate -output=score_html -outfile=kalign.tcs.html
    t_coffee -infile $mafft -evaluate -output=score_html -outfile=mafft.tcs.html
    t_coffee -infile $muscle -evaluate -output=score_html -outfile=muscle.tcs.html
    t_coffee -infile $t_coffee -evaluate -output=score_html -outfile=t_coffee.tcs.html
    t_coffee -infile $msasaBlosum -evaluate -output=score_html -outfile=msasaBlosum.tcs.html
    t_coffee -infile $msasaGlobalMs -evaluate -output=score_html -outfile=msasaGlobalMs.tcs.html
    t_coffee -infile $msasaMatching -evaluate -output=score_html -outfile=msasaMatching.tcs.html
    """
}

workflow {
    runClustalOmega(params.data)
    runKalign(params.data)
    runMAFFT(params.data)
    runMuscle(params.data)
    runTCoffee(params.data)
    runMSASA(params.data)

    computeCoreIndex(
        runClustalOmega.out,
        runKalign.out,
        runMAFFT.out,
        runMuscle.out,
        runTCoffee.out,
        runMSASA.out
    )

    computeTransitiveConsistencyScore(
        runClustalOmega.out,
        runKalign.out,
        runMAFFT.out,
        runMuscle.out,
        runTCoffee.out,
        runMSASA.out
    )
}

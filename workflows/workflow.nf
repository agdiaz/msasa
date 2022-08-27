#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.sequences = '/home/adrian/workspace/msasa/test/small_long.fa'
params.reference = '/home/adrian/workspace/msasa/test/small_short.clustal.fa'
params.outputDir = 'small_long'
params.executions = 5

process runClustalOmega {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path inputFile

    output:
        path "*.clustalo.afa", emit: clustalo_align
        path "*.clustalo.time", emit: clustalo_time

    shell:
    '''
    #!/bin/bash
    RANDOM_MSA_INDEX=$(( ( RANDOM % !{params.executions} ) ))
    echo RANDOM_MSA_INDEX = $RANDOM_MSA_INDEX

    for ((i=0;i<!{params.executions};i++));
    do
        OUT_FILE=result_$i.tmp
        echo OUT_FILE=$OUT_FILE
        /usr/bin/time -p -o $OUT_FILE.time /software/clustalo --force -i !{inputFile} -o $OUT_FILE.afa
    done

    cp result_$RANDOM_MSA_INDEX.tmp.afa result.clustalo.afa
    cp result_$RANDOM_MSA_INDEX.tmp.time result.clustalo.time

    echo BEST=result_$RANDOM_MSA_INDEX.tmp.afa, result_$RANDOM_MSA_INDEX.tmp.time
    '''
}

process runKalign {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path inputFile

    output:
        path "*.kalign.afa", emit: kalign_align
        path "*.kalign.time", emit: kalign_time

    shell:
    '''
    #!/bin/bash
    RANDOM_MSA_INDEX=$(( ( RANDOM % !{params.executions} ) ))
    echo RANDOM_MSA_INDEX = $RANDOM_MSA_INDEX

    for ((i=0;i<!{params.executions};i++));
    do
        OUT_FILE=result_$i.tmp
        echo OUT_FILE=$OUT_FILE

        /usr/bin/time -p -o $OUT_FILE.time /software/kalign -input !{inputFile} -output $OUT_FILE.afa -format fasta -quiet
    done

    cp result_$RANDOM_MSA_INDEX.tmp.afa result.kalign.afa
    cp result_$RANDOM_MSA_INDEX.tmp.time result.kalign.time
    '''
}

process runMAFFT {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path inputFile

    output:
        path "*.mafft.afa", emit: mafft_align
        path "*.mafft.time", emit: mafft_time

    shell:
    '''
    #!/bin/bash
    RANDOM_MSA_INDEX=$(( ( RANDOM % !{params.executions} ) ))
    echo RANDOM_MSA_INDEX = $RANDOM_MSA_INDEX

    for ((i=0;i<!{params.executions};i++));
    do
        OUT_FILE=result_$i.tmp
        echo OUT_FILE=$OUT_FILE

        /usr/bin/time -p -o $OUT_FILE.time /usr/local/bin/mafft --auto !{inputFile} > $OUT_FILE.afa
    done

    cp result_$RANDOM_MSA_INDEX.tmp.afa result.mafft.afa
    cp result_$RANDOM_MSA_INDEX.tmp.time result.mafft.time
    '''
}

process runMuscle {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path inputFile

    output:
        path "*.muscle.afa", emit: muscle_align
        path "*.muscle.time", emit: muscle_time

    shell:
    '''
    #!/bin/bash
    RANDOM_MSA_INDEX=$(( ( RANDOM % !{params.executions} ) ))
    echo RANDOM_MSA_INDEX = $RANDOM_MSA_INDEX

    time --output=example.times ls /
    cat example.times

    for ((i=0;i<!{params.executions};i++));
    do
        OUT_FILE=result_$i.tmp
        echo OUT_FILE=$OUT_FILE

        /usr/bin/time -p -o $OUT_FILE.time /software/muscle -align !{inputFile} -output $OUT_FILE.afa
    done

    cp result_$RANDOM_MSA_INDEX.tmp.afa result.muscle.long_afa
    cp result_$RANDOM_MSA_INDEX.tmp.time result.muscle.time

    fold -w 60 -s result_$RANDOM_MSA_INDEX.tmp.afa > result.muscle.afa
    '''
}

process runTCoffee {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path inputFile

    output:
        path "*.t_coffee.afa", emit: t_coffee_align
        path "*.t_coffee.time", emit: t_coffee_time

    shell:
    '''
    #!/bin/bash

    RANDOM_MSA_INDEX=$(( ( RANDOM % !{params.executions} ) ))
    echo RANDOM_MSA_INDEX = $RANDOM_MSA_INDEX

    for ((i=0;i<!{params.executions};i++));
    do
        OUT_FILE=result_$i.tmp
        echo OUT_FILE=$OUT_FILE

        t_coffee -seq !{inputFile}
        /usr/bin/time -p --output=$OUT_FILE.time t_coffee -other_pg seq_reformat -in !{inputFile.baseName}.aln -output fasta_aln > $OUT_FILE.afa
    done

    cp result_$RANDOM_MSA_INDEX.tmp.afa result.t_coffee.afa
    cp result_$RANDOM_MSA_INDEX.tmp.time result.t_coffee.time

    '''
}

process runMSASA {
    debug true
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path inputFile

    output:
        path "*.msasa_single_ms.afa", emit: msasa_single_ms_align
        path "*.msasa_single_ms.time", emit: msasa_single_ms_time

    shell:
    '''
    LOG_SINGLE_MS_FILENAME=!{inputFile.baseName}_single_ms.log
    LOG_SINGLE_MS_FILENAME_SORTED=$LOG_SINGLE_MS_FILENAME.sorted

    for ((i=0;i<!{params.executions - 3};i++));
    do
        OUTPUT_SINGLE_MS_FILENAME=$i.single_ms.afa
        TIME_SINGLE_MS_FILENAME=$i.single_ms.time

        echo "MSASA Execution # $i started"

        /usr/local/bin/python /usr/src/app/src/msa.py --input !{inputFile} \
            --output $OUTPUT_SINGLE_MS_FILENAME \
            --comparer single_ms \
            --n-iterations 500 \
            --temperature 10 \
            --execution-id $i \
            --engine numpy \
            --optimization min >> $LOG_SINGLE_MS_FILENAME
    done
    echo ORIGINAL
    cat $LOG_SINGLE_MS_FILENAME
    cat $LOG_SINGLE_MS_FILENAME | sort -k 3n --field-separator=";" > $LOG_SINGLE_MS_FILENAME_SORTED

    echo SORTED
    cat $LOG_SINGLE_MS_FILENAME_SORTED

    BEST_SINGLE_MS=$(cat $LOG_SINGLE_MS_FILENAME_SORTED | head -n 1 | cut -d ';' -f1 )
    BEST_SINGLE_MS_ID=$(cat $LOG_SINGLE_MS_FILENAME_SORTED | head -n 1 | cut -d ';' -f4 )
    BEST_SINGLE_MS_TIME=$(cat $LOG_SINGLE_MS_FILENAME_SORTED | head -n 1 | cut -d ';' -f5 )

    echo BEST_SINGLE_MS=$BEST_SINGLE_MS
    echo BEST_SINGLE_MS_TIME=$BEST_SINGLE_MS_TIME

    fold -w 60 -s $BEST_SINGLE_MS > result.msasa_single_ms.afa
    echo $BEST_SINGLE_MS_TIME > result.msasa_single_ms.time
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
        path msasaSingleMs

    output:
        path "*.core-index.html"

    shell:
    """
    t_coffee -infile=$clustalo -output=html -score -outfile=clustalo.core-index.html
    t_coffee -infile=$kalign -output=html -score -outfile=kalign.core-index.html
    t_coffee -infile=$mafft -output=html -score -outfile=mafft.core-index.html
    t_coffee -infile=$muscle -output=html -score -outfile=muscle.core-index.html
    t_coffee -infile=$t_coffee -output=html -score -outfile=t_coffee.core-index.html
    t_coffee -infile=$msasaSingleMs -output=html -score -outfile=msasaSingleMs.core-index.html
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
        path msasaSingleMs

    output:
        path "*.tcs.html"

    shell:
    """
    t_coffee -infile $clustalo -evaluate -output=score_html -outfile=clustalo.tcs.html
    t_coffee -infile $kalign -evaluate -output=score_html -outfile=kalign.tcs.html
    t_coffee -infile $mafft -evaluate -output=score_html -outfile=mafft.tcs.html
    t_coffee -infile $muscle -evaluate -output=score_html -outfile=muscle.tcs.html
    t_coffee -infile $t_coffee -evaluate -output=score_html -outfile=t_coffee.tcs.html
    t_coffee -infile $msasaSingleMs -evaluate -output=score_html -outfile=msasa_single_ms.tcs.html
    """
}

process computeMatrixScore {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path clustalo
        path kalign
        path mafft
        path muscle
        path t_coffee
        path msasaSingleMs

    output:
        path "*.tsv"

    shell:
    """
    CIAlign --infile $clustalo --outfile_stem clustalo --make_similarity_matrix_input --make_similarity_matrix_output
    CIAlign --infile $kalign --outfile_stem kalign --make_similarity_matrix_input --make_similarity_matrix_output
    CIAlign --infile $mafft --outfile_stem mafft --make_similarity_matrix_input --make_similarity_matrix_output
    CIAlign --infile $muscle --outfile_stem muscle --make_similarity_matrix_input --make_similarity_matrix_output
    CIAlign --infile $t_coffee --outfile_stem t_coffee --make_similarity_matrix_input --make_similarity_matrix_output
    CIAlign --infile $msasaSingleMs --outfile_stem msasa_single_ms --make_similarity_matrix_input --make_similarity_matrix_output
    """
}

process computeMumsaOverlapScore {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path clustalo
        path kalign
        path mafft
        path muscle
        path t_coffee
        path msasaSingleMs

    output:
        path "*.overlap"

    shell:
    """
    /software/mumsa-1.0/mumsa -r $clustalo $clustalo $kalign $mafft $t_coffee $msasaSingleMs > results.overlap
    """
}

workflow {
    runClustalOmega(params.sequences)
    runKalign(params.sequences)
    runMAFFT(params.sequences)
    runMuscle(params.sequences)
    runTCoffee(params.sequences)
    runMSASA(params.sequences)

    computeCoreIndex(
        runClustalOmega.out.clustalo_align,
        runKalign.out.kalign_align,
        runMAFFT.out.mafft_align,
        runMuscle.out.muscle_align,
        runTCoffee.out.t_coffee_align,
        runMSASA.out.msasa_single_ms_align
    )

    computeTransitiveConsistencyScore(
        runClustalOmega.out.clustalo_align,
        runKalign.out.kalign_align,
        runMAFFT.out.mafft_align,
        runMuscle.out.muscle_align,
        runTCoffee.out.t_coffee_align,
        runMSASA.out.msasa_single_ms_align
    )

    computeMatrixScore(
        runClustalOmega.out.clustalo_align,
        runKalign.out.kalign_align,
        runMAFFT.out.mafft_align,
        runMuscle.out.muscle_align,
        runTCoffee.out.t_coffee_align,
        runMSASA.out.msasa_single_ms_align
    )

    computeMumsaOverlapScore(
        runClustalOmega.out.clustalo_align,
        runKalign.out.kalign_align,
        runMAFFT.out.mafft_align,
        runMuscle.out.muscle_align,
        runTCoffee.out.t_coffee_align,
        runMSASA.out.msasa_single_ms_align
    )
}

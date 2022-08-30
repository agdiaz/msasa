#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.sequences = '/home/adrian/workspace/msasa/database/bb3_release/RV20/BB20001.tfa'
params.reference = '/home/adrian/workspace/msasa/database/bb3_release/RV20/BB20001'
params.outputDir = 'BB20001'

params.executions = 3

process runClustalOmega {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path inputFile

    output:
        path "clustalo.afa", emit: clustalo_align
        path "clustalo.time", emit: clustalo_time

    shell:
    '''
    #!/bin/bash
    RANDOM_MSA_INDEX=$(( ( RANDOM % !{params.executions} ) ))

    for ((i=0;i<!{params.executions};i++));
    do
        OUT_FILE=result_$i
        echo OUT_FILE=$OUT_FILE
        /usr/bin/time -p -o $OUT_FILE.time /software/clustalo --force -i !{inputFile} -o $OUT_FILE.afa
    done

    cp result_$RANDOM_MSA_INDEX.afa clustalo.afa
    cp result_$RANDOM_MSA_INDEX.time clustalo.time
    '''
}

process runKalign {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path inputFile

    output:
        path "kalign.afa", emit: kalign_align
        path "kalign.time", emit: kalign_time

    shell:
    '''
    #!/bin/bash
    RANDOM_MSA_INDEX=$(( ( RANDOM % !{params.executions} ) ))

    for ((i=0;i<!{params.executions};i++));
    do
        OUT_FILE=result_$i
        echo OUT_FILE=$OUT_FILE

        /usr/bin/time -p -o $OUT_FILE.time /software/kalign -input !{inputFile} -output $OUT_FILE.afa -format fasta -quiet
    done

    cp result_$RANDOM_MSA_INDEX.afa kalign.afa
    cp result_$RANDOM_MSA_INDEX.time kalign.time
    '''
}

process runMAFFT {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path inputFile

    output:
        path "mafft.afa", emit: mafft_align
        path "mafft.time", emit: mafft_time

    shell:
    '''
    #!/bin/bash
    RANDOM_MSA_INDEX=$(( ( RANDOM % !{params.executions} ) ))

    for ((i=0;i<!{params.executions};i++));
    do
        OUT_FILE=result_$i
        echo OUT_FILE=$OUT_FILE

        /usr/bin/time -p -o $OUT_FILE.time /usr/local/bin/mafft --auto !{inputFile} > $OUT_FILE.afa
    done

    cp result_$RANDOM_MSA_INDEX.afa mafft.afa
    cp result_$RANDOM_MSA_INDEX.time mafft.time
    '''
}

process runMuscle {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path inputFile

    output:
        path "muscle.afa", emit: muscle_align
        path "muscle.time", emit: muscle_time

    shell:
    '''
    #!/bin/bash
    RANDOM_MSA_INDEX=$(( ( RANDOM % !{params.executions} ) ))

    for ((i=0;i<!{params.executions};i++));
    do
        OUT_FILE=result_$i
        echo OUT_FILE=$OUT_FILE

        /usr/bin/time -p -o $OUT_FILE.time /software/muscle -align !{inputFile} -output $OUT_FILE.afa
    done

    fold -w 60 -s result_$RANDOM_MSA_INDEX.afa > muscle.afa
    cp result_$RANDOM_MSA_INDEX.time muscle.time
    '''
}

process runTCoffee {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path inputFile

    output:
        path "t_coffee.afa", emit: t_coffee_align
        path "t_coffee.time", emit: t_coffee_time

    shell:
    '''
    #!/bin/bash

    RANDOM_MSA_INDEX=$(( ( RANDOM % !{params.executions} ) ))
    echo RANDOM_MSA_INDEX = $RANDOM_MSA_INDEX

    for ((i=0;i<!{params.executions};i++));
    do
        OUT_FILE=result_$i
        echo OUT_FILE=$OUT_FILE

        t_coffee -seq !{inputFile}
        /usr/bin/time -p --output=$OUT_FILE.time t_coffee -other_pg seq_reformat -in !{inputFile.baseName}.aln -output fasta_aln > $OUT_FILE.afa
    done

    cp result_$RANDOM_MSA_INDEX.afa t_coffee.afa
    cp result_$RANDOM_MSA_INDEX.time t_coffee.time
    '''
}

process runMSASA {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'
    cpus 12
    errorStrategy 'terminate'

    input:
        path inputFile

    output:
        path "msasa_single_ms.afa", emit: msasa_single_ms_align
        path "msasa_single_ms.time", emit: msasa_single_ms_time

    shell:
    '''
    LOG_SINGLE_MS_FILENAME=!{inputFile.baseName}_single_ms.log
    LOG_SINGLE_MS_FILENAME_SORTED=$LOG_SINGLE_MS_FILENAME.sorted

    for ((i=0;i<!{params.executions};i++));
    do
        OUTPUT_SINGLE_MS_FILENAME=$i.single_ms.afa
        TIME_SINGLE_MS_FILENAME=$i.single_ms.time

        /usr/local/bin/python /usr/src/app/src/msa.py --input !{inputFile} \
            --output $OUTPUT_SINGLE_MS_FILENAME \
            --comparer single_ms \
            --n-iterations 5000 \
            --temperature 50 \
            --execution-id $i \
            --engine numpy \
            --optimization min >> $LOG_SINGLE_MS_FILENAME
    done

    cat $LOG_SINGLE_MS_FILENAME | sort -k 3n --field-separator=";" > $LOG_SINGLE_MS_FILENAME_SORTED

    BEST_SINGLE_MS=$(cat $LOG_SINGLE_MS_FILENAME_SORTED | head -n 1 | cut -d ';' -f1 )
    BEST_SINGLE_MS_ID=$(cat $LOG_SINGLE_MS_FILENAME_SORTED | head -n 1 | cut -d ';' -f4 )
    BEST_SINGLE_MS_TIME=$(cat $LOG_SINGLE_MS_FILENAME_SORTED | head -n 1 | cut -d ';' -f5 )

    fold -w 60 -s $BEST_SINGLE_MS > msasa_single_ms.afa
    echo $BEST_SINGLE_MS_TIME > msasa_single_ms.time
    '''
}

process runMSASABlosum {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'
    cpus 12
    errorStrategy 'terminate'

    input:
        path inputFile

    output:
        path "msasa_blosum.afa", emit: msasa_blosum_align
        path "msasa_blosum.time", emit: msasa_blosum_time

    shell:
    '''
    LOG_SINGLE_BLOSUM_FILENAME=!{inputFile.baseName}_single_blosum.log
    LOG_SINGLE_BLOSUM_FILENAME_SORTED=$LOG_SINGLE_BLOSUM_FILENAME.sorted

    for ((i=0;i<!{params.executions};i++));
    do
        OUTPUT_SINGLE_BLOSUM_FILENAME=$i.single_blosum.afa
        TIME_SINGLE_BLOSUM_FILENAME=$i.single_blosum.time

        /usr/local/bin/python /usr/src/app/src/msa.py --input !{inputFile} \
            --output $OUTPUT_SINGLE_BLOSUM_FILENAME \
            --comparer single_blosum \
            --n-iterations 5000 \
            --temperature 50 \
            --execution-id $i \
            --engine numpy \
            --optimization min >> $LOG_SINGLE_BLOSUM_FILENAME
    done
    cat $LOG_SINGLE_BLOSUM_FILENAME | sort -k 3n --field-separator=";" > $LOG_SINGLE_BLOSUM_FILENAME_SORTED

    BEST_SINGLE_BLOSUM=$(cat $LOG_SINGLE_BLOSUM_FILENAME_SORTED | head -n 1 | cut -d ';' -f1 )
    BEST_SINGLE_BLOSUM_ID=$(cat $LOG_SINGLE_BLOSUM_FILENAME_SORTED | head -n 1 | cut -d ';' -f4 )
    BEST_SINGLE_BLOSUM_TIME=$(cat $LOG_SINGLE_BLOSUM_FILENAME_SORTED | head -n 1 | cut -d ';' -f5 )

    fold -w 60 -s $BEST_SINGLE_BLOSUM > msasa_blosum.afa
    echo $BEST_SINGLE_BLOSUM_TIME > msasa_blosum.time
    '''
}

process runMSASAMatching {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'
    cpus 12
    errorStrategy 'terminate'

    input:
        path inputFile

    output:
        path "msasa_single_matching.afa", emit: msasa_matching_align
        path "msasa_single_matching.time", emit: msasa_matching_time

    shell:
    '''
    LOG_SINGLE_MATCHING_FILENAME=!{inputFile.baseName}_single_matching.log
    LOG_SINGLE_MATCHING_FILENAME_SORTED=$LOG_SINGLE_MATCHING_FILENAME.sorted

    for ((i=0;i<!{params.executions};i++));
    do
        OUTPUT_SINGLE_MATCHING_FILENAME=$i.single_matching.afa
        TIME_SINGLE_MATCHING_FILENAME=$i.single_matching.time

        /usr/local/bin/python /usr/src/app/src/msa.py --input !{inputFile} \
            --output $OUTPUT_SINGLE_MATCHING_FILENAME \
            --comparer single_matching \
            --n-iterations 5000 \
            --temperature 50 \
            --execution-id $i \
            --engine numpy \
            --optimization min >> $LOG_SINGLE_MATCHING_FILENAME
    done
    cat $LOG_SINGLE_MATCHING_FILENAME | sort -k 3n --field-separator=";" > $LOG_SINGLE_MATCHING_FILENAME_SORTED

    BEST_SINGLE_MATCHING=$(cat $LOG_SINGLE_MATCHING_FILENAME_SORTED | head -n 1 | cut -d ';' -f1 )
    BEST_SINGLE_MATCHING_ID=$(cat $LOG_SINGLE_MATCHING_FILENAME_SORTED | head -n 1 | cut -d ';' -f4 )
    BEST_SINGLE_MATCHING_TIME=$(cat $LOG_SINGLE_MATCHING_FILENAME_SORTED | head -n 1 | cut -d ';' -f5 )

    fold -w 60 -s $BEST_SINGLE_MATCHING > msasa_single_matching.afa
    echo $BEST_SINGLE_MATCHING_TIME > msasa_single_matching.time
    '''
}

process convertFastaToMSF {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path clustalo
        path kalign
        path mafft
        path muscle
        path t_coffee
        path msasaSingleMs
        path msasaSingleBlosum
        path msasaMatchingAlign

    output:
        path "clustalo.msf", emit: clustaloMsf
        path "kalign.msf", emit: kalignMsf
        path "mafft.msf", emit: mafftMsf
        path "muscle.msf", emit: muscleMsf
        path "t_coffee.msf", emit: t_coffeeMsf
        path "msasaSingleMs.msf", emit: msasaSingleMsMsf
        path "msasaSingleBlosum.msf", emit: msasaSingleBlosumMsf
        path "msasaMatchingAlign.msf", emit: msasaMatchingAlignMsf

    shell:
    '''
    seqlim -outfmt msf -o clustalo.msf cnvt !{clustalo}
    seqlim -outfmt msf -o kalign.msf cnvt !{kalign}
    seqlim -outfmt msf -o mafft.msf cnvt !{mafft}
    seqlim -outfmt msf -o muscle.msf cnvt !{muscle}
    seqlim -outfmt msf -o t_coffee.msf cnvt !{t_coffee}
    seqlim -outfmt msf -o msasaSingleMs.msf cnvt !{msasaSingleMs}
    seqlim -outfmt msf -o msasaSingleBlosum.msf cnvt !{msasaSingleBlosum}
    seqlim -outfmt msf -o msasaMatchingAlign.msf cnvt !{msasaMatchingAlign}
    '''
}

process convertReferenceFromMsfToFasta {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    output:
        path "reference.fasta", emit: reference_fasta

    shell:
    """
    seqlim cnvt ${params.reference}.msf -infmt msf -outfmt fasta -o reference.fasta
    """
}

process computeCoreIndexClustalo {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.core-index.html"

    shell:
    """
    t_coffee -infile=$alignment -output=html -score -outfile=${alignment.baseName}.core-index.html
    """
}

process computeCoreIndexKalign {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.core-index.html"

    shell:
    """
    t_coffee -infile=$alignment -output=html -score -outfile=${alignment.baseName}.core-index.html
    """
}

process computeCoreIndexMafft {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.core-index.html"

    shell:
    """
    t_coffee -infile=$alignment -output=html -score -outfile=${alignment.baseName}.core-index.html
    """
}

process computeCoreIndexMuscle {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.core-index.html"

    shell:
    """
    t_coffee -infile=$alignment -output=html -score -outfile=${alignment.baseName}.core-index.html
    """
}

process computeCoreIndexTCoffee {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.core-index.html"

    shell:
    """
    t_coffee -infile=$alignment -output=html -score -outfile=${alignment.baseName}.core-index.html
    """
}

process computeCoreIndexMsasa {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.core-index.html"

    shell:
    """
    t_coffee -infile=$alignment -output=html -score -outfile=${alignment.baseName}.core-index.html
    """
}

process computeCoreIndexMsasaBlosum {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.core-index.html"

    shell:
    """
    t_coffee -infile=$alignment -output=html -score -outfile=${alignment.baseName}.core-index.html
    """
}

process computeCoreIndexMsasaMatching {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.core-index.html"

    shell:
    """
    t_coffee -infile=$alignment -output=html -score -outfile=${alignment.baseName}.core-index.html
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
        path msasaSingleBlosum
        path msasaMatchingAlign

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
    t_coffee -infile $msasaSingleBlosum -evaluate -output=score_html -outfile=msasa_single_blosum.tcs.html
    t_coffee -infile $msasaMatchingAlign -evaluate -output=score_html -outfile=msasa_single_matching.tcs.html
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
        path msasaSingleBlosum
        path msasaMatchingAlign

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
    CIAlign --infile $msasaSingleBlosum --outfile_stem msasa_single_blosum --make_similarity_matrix_input --make_similarity_matrix_output
    CIAlign --infile $msasaMatchingAlign --outfile_stem msasa_single_matching --make_similarity_matrix_input --make_similarity_matrix_output
    """
}

process computeMumsaOverlapScore {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path reference_fasta
        path clustalo
        path kalign
        path mafft
        path muscle
        path t_coffee
        path msasaSingleMs
        path msasaSingleBlosum
        path msasaMatchingAlign

    output:
        path "results.overlap.txt"

    shell:
    """
    /software/mumsa-1.0/mumsa -r -q $reference_fasta $clustalo $kalign $mafft $t_coffee $msasaSingleMs $msasaSingleBlosum $msasaMatchingAlign > results.overlap.txt
    """
}

process computeBaliScoreClustal {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.baliscore.txt"

    shell:
    '''
    /software/bali-score/target/release/bali-score -t !{alignment} -r !{params.reference}.xml -o !{alignment}.baliscore.txt
    '''
}

process computeBaliScoreKalign {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.baliscore.txt"

    shell:
    '''
    /software/bali-score/target/release/bali-score -t !{alignment} -r !{params.reference}.xml -o !{alignment}.baliscore.txt
    '''
}

process computeBaliScoreMafft {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.baliscore.txt"

    shell:
    '''
    /software/bali-score/target/release/bali-score -t !{alignment} -r !{params.reference}.xml -o !{alignment}.baliscore.txt
    '''
}

process computeBaliScoreMuscle {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.baliscore.txt"

    shell:
    '''
    /software/bali-score/target/release/bali-score -t !{alignment} -r !{params.reference}.xml -o !{alignment}.baliscore.txt
    '''
}

process computeBaliScoreTCoffee {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.baliscore.txt"

    shell:
    '''
    /software/bali-score/target/release/bali-score -t !{alignment} -r !{params.reference}.xml -o !{alignment}.baliscore.txt
    '''
}

process computeBaliScoreMsasa {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.baliscore.txt"

    shell:
    '''
    /software/bali-score/target/release/bali-score -t !{alignment} -r !{params.reference}.xml -o !{alignment}.baliscore.txt
    '''
}


process computeBaliScoreMsasaBlosum {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.baliscore.txt"

    shell:
    '''
    /software/bali-score/target/release/bali-score -t !{alignment} -r !{params.reference}.xml -o !{alignment}.baliscore.txt
    '''
}

process computeBaliScoreMsasaMatching {
    publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

    input:
        path alignment

    output:
        path "*.baliscore.txt"

    shell:
    '''
    /software/bali-score/target/release/bali-score -t !{alignment} -r !{params.reference}.xml -o !{alignment}.baliscore.txt
    '''
}


workflow {
    convertReferenceFromMsfToFasta()

    runClustalOmega(params.sequences)
    runKalign(params.sequences)
    runMAFFT(params.sequences)
    runMuscle(params.sequences)
    runTCoffee(params.sequences)
    runMSASA(params.sequences)
    runMSASABlosum(params.sequences)
    runMSASAMatching(params.sequences)

    computeCoreIndexClustalo(runClustalOmega.out.clustalo_align)
    computeCoreIndexKalign(runKalign.out.kalign_align)
    computeCoreIndexMafft(runMAFFT.out.mafft_align)
    computeCoreIndexMuscle(runMuscle.out.muscle_align)
    computeCoreIndexTCoffee(runTCoffee.out.t_coffee_align)
    computeCoreIndexMsasa(runMSASA.out.msasa_single_ms_align)
    computeCoreIndexMsasaBlosum(runMSASABlosum.out.msasa_blosum_align)
    computeCoreIndexMsasaMatching(runMSASAMatching.out.msasa_matching_align)

    // computeTransitiveConsistencyScore(
    //     runClustalOmega.out.clustalo_align,
    //     runKalign.out.kalign_align,
    //     runMAFFT.out.mafft_align,
    //     runMuscle.out.muscle_align,
    //     runTCoffee.out.t_coffee_align,
    //     runMSASA.out.msasa_single_ms_align,
    //     runMSASABlosum.out.msasa_blosum_align,
    //     runMSASAMatching.out.msasa_matching_align
    // )

    // computeMatrixScore(
    //     runClustalOmega.out.clustalo_align,
    //     runKalign.out.kalign_align,
    //     runMAFFT.out.mafft_align,
    //     runMuscle.out.muscle_align,
    //     runTCoffee.out.t_coffee_align,
    //     runMSASA.out.msasa_single_ms_align,
    //     runMSASABlosum.out.msasa_blosum_align,
    //     runMSASAMatching.out.msasa_matching_align
    // )

    computeMumsaOverlapScore(
        convertReferenceFromMsfToFasta.out.reference_fasta,
        runClustalOmega.out.clustalo_align,
        runKalign.out.kalign_align,
        runMAFFT.out.mafft_align,
        runMuscle.out.muscle_align,
        runTCoffee.out.t_coffee_align,
        runMSASA.out.msasa_single_ms_align,
        runMSASABlosum.out.msasa_blosum_align,
        runMSASAMatching.out.msasa_matching_align
    )

    computeBaliScoreClustal(runClustalOmega.out.clustalo_align)
    computeBaliScoreKalign(runKalign.out.kalign_align)
    computeBaliScoreMafft(runMAFFT.out.mafft_align)
    // computeBaliScoreMuscle(runMuscle.out.muscle_align)
    computeBaliScoreTCoffee(runTCoffee.out.t_coffee_align)
    computeBaliScoreMsasa(runMSASA.out.msasa_single_ms_align)
    computeBaliScoreMsasaBlosum(runMSASABlosum.out.msasa_blosum_align)
    computeBaliScoreMsasaMatching(runMSASAMatching.out.msasa_matching_align)
}

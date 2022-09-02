#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.sequences = '/home/adrian/workspace/msasa/database/bb3_release/RV12/BB12001.tfa'
params.reference = '/home/adrian/workspace/msasa/database/bb3_release/RV12/BB12001.msf'
params.referenceXml = '/home/adrian/workspace/msasa/database/bb3_release/RV12/BB12001.xml'
params.executions = 30

process runClustalOmega {
    publishDir "$baseDir/results/$sequences.simpleName", mode: 'copy'

    input:
        tuple val(executionId), path(sequences)

    output:
        tuple path("$executionId.${sequences.simpleName}.clustalo.afa"), path("$executionId.${sequences.simpleName}.clustalo.time")

    shell:
    '''
    #!/bin/bash
    /usr/bin/time -p -o !{executionId}.!{sequences.simpleName}.clustalo.time /software/clustalo --force -i !{sequences} -o !{executionId}.!{sequences.simpleName}.clustalo.afa
    '''
}

process runKalign {
    publishDir "$baseDir/results/$sequences.simpleName", mode: 'copy'

    input:
        tuple val(executionId), path(sequences)

    output:
        tuple path("$executionId.${sequences.simpleName}.kalign.afa"), path("$executionId.${sequences.simpleName}.kalign.time")

    shell:
    '''
    #!/bin/bash
    /usr/bin/time -p -o !{executionId}.!{sequences.simpleName}.kalign.time /software/kalign -input !{sequences} -output !{executionId}.!{sequences.simpleName}.kalign.afa -format fasta -quiet
    '''
}

process runMAFFT {
    publishDir "$baseDir/results/$sequences.simpleName", mode: 'copy'

    input:
        tuple val(executionId), path(sequences)

    output:
        tuple path("$executionId.${sequences.simpleName}.mafft.afa"), path("$executionId.${sequences.simpleName}.mafft.time")

    shell:
    '''
    #!/bin/bash
    /usr/bin/time -p -o !{executionId}.!{sequences.simpleName}.mafft.time /usr/local/bin/mafft --auto !{sequences} > !{executionId}.!{sequences.simpleName}.mafft.afa
    '''
}

process runMuscle {
    publishDir "$baseDir/results/$sequences.simpleName", mode: 'copy'

    input:
        tuple val(executionId), path(sequences)

    output:
        tuple path("$executionId.${sequences.simpleName}.muscle.afa"), path("$executionId.${sequences.simpleName}.muscle.time")

    shell:
    '''
    /usr/bin/time -p -o !{executionId}.!{sequences.simpleName}.muscle.time /software/muscle -align !{sequences} -output temp.afa

    fold -w 60 -s temp.afa > !{executionId}.!{sequences.simpleName}.muscle.afa
    '''
}

process runTCoffee {
    publishDir "$baseDir/results/$sequences.simpleName", mode: 'copy'

    input:
        tuple val(executionId), path(sequences)

    output:
        tuple path("$executionId.${sequences.simpleName}.t_coffee.afa"), path("$executionId.${sequences.simpleName}.t_coffee.time")

    shell:
    '''
    #!/bin/bash
    t_coffee -seq !{sequences}
    /usr/bin/time -p --output=!{executionId}.!{sequences.simpleName}.t_coffee.time t_coffee -other_pg seq_reformat -in !{sequences.simpleName}.aln -output fasta_aln > !{executionId}.!{sequences.simpleName}.t_coffee.afa
    '''
}

process runMSASA {
    publishDir "$baseDir/results/$sequences.simpleName", mode: 'copy'
    cpus 12

    input:
        tuple val(executionId), path(sequences)

    output:
        tuple path("$executionId.${sequences.simpleName}.single_ms.afa"), path("$executionId.${sequences.simpleName}.single_ms.time")

    shell:
    '''
    #!/bin/bash
    /usr/local/bin/python /usr/src/app/src/msa.py --input !{sequences} \
        --output !{executionId}.!{sequences.simpleName}.single_ms.afa \
        --comparer single_ms \
        --n-iterations 100 \
        --temperature 10 \
        --execution-id 1 > !{executionId}.!{sequences.simpleName}.single_ms.time
    '''
}

process runMSASABlosum {
    publishDir "$baseDir/results/$sequences.simpleName", mode: 'copy'
    cpus 12

    input:
        tuple val(executionId), path(sequences)

    output:
        tuple path("$executionId.${sequences.simpleName}.single_blosum.afa"), path("$executionId.${sequences.simpleName}.single_blosum.time")

    shell:
    '''
    #!/bin/bash
    /usr/local/bin/python /usr/src/app/src/msa.py --input !{sequences} \
        --output !{executionId}.!{sequences.simpleName}.single_blosum.afa \
        --comparer single_blosum \
        --n-iterations 100 \
        --temperature 20 \
        --execution-id 1 > !{executionId}.!{sequences.simpleName}.single_blosum.time
    '''
}

// process convertFastaToMSF {
//     publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

//     input:
//         path clustalo
//         path kalign
//         path mafft
//         path muscle
//         path t_coffee
//         path msasaSingleMs
//         path msasaSingleBlosum

//     output:
//         path "clustalo.msf", emit: clustaloMsf
//         path "kalign.msf", emit: kalignMsf
//         path "mafft.msf", emit: mafftMsf
//         path "muscle.msf", emit: muscleMsf
//         path "t_coffee.msf", emit: t_coffeeMsf
//         path "msasaSingleMs.msf", emit: msasaSingleMsMsf
//         path "msasaSingleBlosum.msf", emit: msasaSingleBlosumMsf

//     shell:
//     '''
//     seqlim -outfmt msf -o clustalo.msf cnvt !{clustalo}
//     seqlim -outfmt msf -o kalign.msf cnvt !{kalign}
//     seqlim -outfmt msf -o mafft.msf cnvt !{mafft}
//     seqlim -outfmt msf -o muscle.msf cnvt !{muscle}
//     seqlim -outfmt msf -o t_coffee.msf cnvt !{t_coffee}
//     seqlim -outfmt msf -o msasaSingleMs.msf cnvt !{msasaSingleMs}
//     seqlim -outfmt msf -o msasaSingleBlosum.msf cnvt !{msasaSingleBlosum}
//     '''
// }

process convertReferenceFromMsfToFasta {
    output:
        path "reference.fasta", emit: fasta

    shell:
    '''
    seqlim cnvt !{params.reference} -infmt msf -outfmt fasta -o reference.fasta
    '''
}

// process computeCoreIndexClustalo {
//     publishDir "$baseDir/results/$alignment.simpleName", mode: 'copy'

//     input:
//         path alignment

//     output:
//         path "*.core-index.html"

//     shell:
//     """
//     t_coffee -infile=$alignment -output=html -score -outfile=${alignment.simpleName}.coreindex.html
//     """
// }

// process computeCoreIndexKalign {
//     publishDir "$baseDir/results/$alignment.simpleName", mode: 'copy'

//     input:
//         path alignment

//     output:
//         path "*.core-index.html"

//     shell:
//     """
//     t_coffee -infile=$alignment -output=html -score -outfile=${alignment.simpleName}.coreindex.html
//     """
// }

// process computeCoreIndexMafft {
//     publishDir "$baseDir/results/$alignment.simpleName", mode: 'copy'

//     input:
//         path alignment

//     output:
//         path "*.core-index.html"

//     shell:
//     """
//     t_coffee -infile=$alignment -output=html -score -outfile=${alignment.simpleName}.coreindex.html
//     """
// }

// process computeCoreIndexMuscle {
//     publishDir "$baseDir/results/$alignment.simpleName", mode: 'copy'

//     input:
//         path alignment

//     output:
//         path "*.core-index.html"

//     shell:
//     """
//     t_coffee -infile=$alignment -output=html -score -outfile=${alignment.simpleName}.coreindex.html
//     """
// }

// process computeCoreIndexTCoffee {
//     publishDir "$baseDir/results/$alignment.simpleName", mode: 'copy'

//     input:
//         path alignment

//     output:
//         path "*.core-index.html"

//     shell:
//     """
//     t_coffee -infile=$alignment -output=html -score -outfile=${alignment.simpleName}.coreindex.html
//     """
// }

// process computeCoreIndexMsasa {
//     publishDir "$baseDir/results/$alignment.simpleName", mode: 'copy'

//     input:
//         path alignment

//     output:
//         path "*.core-index.html"

//     shell:
//     """
//     t_coffee -infile=$alignment -output=html -score -outfile=${alignment.simpleName}.coreindex.html
//     """
// }

// process computeCoreIndexMsasaBlosum {
//     publishDir "$baseDir/results/$alignment.simpleName", mode: 'copy'

//     input:
//         path alignment

//     output:
//         path "*.core-index.html"

//     shell:
//     """
//     t_coffee -infile=$alignment -output=html -score -outfile=${alignment.simpleName}.coreindex.html
//     """
// }


// process computeTransitiveConsistencyScore {
//     publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

//     input:
//         path clustalo
//         path kalign
//         path mafft
//         path muscle
//         path t_coffee
//         path msasaSingleMs
//         path msasaSingleBlosum

//     output:
//         path "*.tcs.html"

//     shell:
//     """
//     t_coffee -infile $clustalo -evaluate -output=score_html -outfile=clustalo.tcs.html
//     t_coffee -infile $kalign -evaluate -output=score_html -outfile=kalign.tcs.html
//     t_coffee -infile $mafft -evaluate -output=score_html -outfile=mafft.tcs.html
//     t_coffee -infile $muscle -evaluate -output=score_html -outfile=muscle.tcs.html
//     t_coffee -infile $t_coffee -evaluate -output=score_html -outfile=t_coffee.tcs.html
//     t_coffee -infile $msasaSingleMs -evaluate -output=score_html -outfile=msasa_single_ms.tcs.html
//     t_coffee -infile $msasaSingleBlosum -evaluate -output=score_html -outfile=msasa_single_blosum.tcs.html
//     """
// }

// process computeMatrixScore {
//     publishDir "$baseDir/results/$params.outputDir", mode: 'copy'

//     input:
//         path clustalo
//         path kalign
//         path mafft
//         path muscle
//         path t_coffee
//         path msasaSingleMs
//         path msasaSingleBlosum

//     output:
//         path "*.tsv"

//     shell:
//     """
//     CIAlign --infile $clustalo --outfile_stem clustalo --make_similarity_matrix_input --make_similarity_matrix_output
//     CIAlign --infile $kalign --outfile_stem kalign --make_similarity_matrix_input --make_similarity_matrix_output
//     CIAlign --infile $mafft --outfile_stem mafft --make_similarity_matrix_input --make_similarity_matrix_output
//     CIAlign --infile $muscle --outfile_stem muscle --make_similarity_matrix_input --make_similarity_matrix_output
//     CIAlign --infile $t_coffee --outfile_stem t_coffee --make_similarity_matrix_input --make_similarity_matrix_output
//     CIAlign --infile $msasaSingleMs --outfile_stem msasa_single_ms --make_similarity_matrix_input --make_similarity_matrix_output
//     CIAlign --infile $msasaSingleBlosum --outfile_stem msasa_single_blosum --make_similarity_matrix_input --make_similarity_matrix_output
//     """
// }

process computeMumsaOverlapScore {
    publishDir "$baseDir/results/$clustalo.simpleName", mode: 'copy'

    input:
        path reference_fasta
        path clustalo
        path kalign
        path mafft
        path muscle
        path t_coffee
        path msasaSingleMs
        path msasaSingleBlosum

    output:
        path "*.overlap.txt"

    shell:
    """
    /software/mumsa-1.0/mumsa -r -q $reference_fasta $clustalo $kalign $mafft $t_coffee $msasaSingleMs $msasaSingleBlosum > ${clustalo.simpleName}.overlap.txt
    """
}

process computeBaliScore {
    publishDir "$baseDir/results/$alignment.simpleName", mode: 'copy'

    input:
        path alignment

    output:
        path "*.baliscore.txt"

    shell:
    '''
    /software/bali-score/target/release/bali-score -t !{alignment} -r !{params.referenceXml} -o !{alignment.baseName}.baliscore.txt
    '''
}

process getClustalSample {

    input:
        tuple path(x), path(y)
    output:
        path(x), emit: sampleAlignment
        path(y), emit: sample_timing

    script:
    """
    echo $x
    echo $y
    """
}

process getKalignSample {
    input:
        tuple path(x), path(y)
    output:
        path(x), emit: sampleAlignment
        path(y), emit: sample_timing

    script:
    """
    echo $x
    echo $y
    """
}

process getMafftSample {
    input:
        tuple path(x), path(y)
    output:
        path(x), emit: sampleAlignment
        path(y), emit: sample_timing

    script:
    """
    echo $x
    echo $y
    """
}

process getMuscleSample {
    input:
        tuple path(x), path(y)
    output:
        path(x), emit: sampleAlignment
        path(y), emit: sample_timing

    script:
    """
    echo $x
    echo $y
    """
}

process getTCoffeeSample {
    input:
        tuple path(x), path(y)
    output:
        path(x), emit: sampleAlignment
        path(y), emit: sample_timing

    script:
    """
    echo $x
    echo $y
    """
}

process getSingleMSSample {
    input:
        tuple path(x), path(y)
    output:
        path(x), emit: sampleAlignment
        path(y), emit: sample_timing

    script:
    """
    echo $x
    echo $y
    """
}

process getSingleBlosumSample {
    input:
        tuple path(x), path(y)
    output:
        path(x), emit: sampleAlignment
        path(y), emit: sample_timing

    script:
    """
    echo $x
    echo $y
    """
}

workflow runMSAtools {
    take: sequences
    main:
        runClustalOmega(sequences)
        runKalign(sequences)
        runMAFFT(sequences)
        runMuscle(sequences)
        runTCoffee(sequences)
        runMSASA(sequences)
        runMSASABlosum(sequences)

    emit:
        clustalAlignment = runClustalOmega.out
        kalignAlignment = runKalign.out
        mafftAlignment = runMAFFT.out
        muscleAlignment = runMuscle.out
        tcoffeeAlignment = runTCoffee.out
        singleMsAlignment = runMSASA.out
        singleBlosumAlignment = runMSASABlosum.out

}

workflow computeAnalytics {
    take:
        clustalAlignment
        kalignAlignment
        mafftAlignment
        muscleAlignment
        tcoffeeAlignment
        singleMsAlignment
        singleBlosumAlignment
    main:
        convertReferenceFromMsfToFasta()

        computeMumsaOverlapScore(
            convertReferenceFromMsfToFasta.out.fasta,
            clustalAlignment,
            kalignAlignment,
            mafftAlignment,
            muscleAlignment,
            tcoffeeAlignment,
            singleMsAlignment,
            singleBlosumAlignment,
        )
}

workflow {
    sequencesChannel = Channel.from( 1 .. params.iterations ).map { it -> [ number: it, sequences: params.sequences ] }

    runClustalOmega(sequencesChannel)
    runKalign(sequencesChannel)
    runMAFFT(sequencesChannel)
    runMuscle(sequencesChannel)
    runTCoffee(sequencesChannel)
    runMSASA(sequencesChannel)
    runMSASABlosum(sequencesChannel)

    getClustalSample(runClustalOmega.out.randomSample(1))
    getKalignSample(runKalign.out.randomSample(1))
    getMafftSample(runMAFFT.out.randomSample(1))
    getMuscleSample(runMuscle.out.randomSample(1))
    getTCoffeeSample(runTCoffee.out.randomSample(1))
    getSingleMSSample(runMSASA.out.randomSample(1))
    getSingleBlosumSample(runMSASABlosum.out.randomSample(1))

    alignments = getClustalSample.out.sampleAlignment.concat(
        getKalignSample.out.sampleAlignment,
        getMafftSample.out.sampleAlignment,
        // getMuscleSample.out.sampleAlignment,
        getTCoffeeSample.out.sampleAlignment,
        getSingleMSSample.out.sampleAlignment,
        getSingleBlosumSample.out.sampleAlignment
    )

    computeAnalytics(
        getClustalSample.out.sampleAlignment,
        getKalignSample.out.sampleAlignment,
        getMafftSample.out.sampleAlignment,
        getMuscleSample.out.sampleAlignment,
        getTCoffeeSample.out.sampleAlignment,
        getSingleMSSample.out.sampleAlignment,
        getSingleBlosumSample.out.sampleAlignment
    )

    computeBaliScore(alignments)
}

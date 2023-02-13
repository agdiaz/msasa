// nextflow run tesis.nf -resume -c tesis.conf --tfaFile ../dataset_testing/BB11004.tfa --referenceXml ../dataset_testing/BB11004.xml --reference ../dataset_testing/BB11004.msf

params.tfaFile = ""
params.reference = ""
params.referenceXml = ""

params.iterations = 5

targetFile = file(params.tfaFile)
referenceFile = file(params.reference)

include {
    predictMuscle;
    predictClustalOmega;
    predictMafft;
    predictTCoffee;
    predictHmmer;
    predictKAlign;
    sortPredictedAlignment;
} from "./modules/predictors/main.nf"

include {
    computeBaliScore as computeBaliScoreMuscle;
    computeBaliScore as computeBaliScoreClustal;
    computeBaliScore as computeBaliScoreMafft;
    computeBaliScore as computeBaliScoreTCoffee;
    computeBaliScore as computeBaliScoreHmmer;
    computeBaliScore as computeBaliScoreKalign;
    computeMumsaOverlapScore;
    exportAlignments;
} from "./modules/analysis/main.nf"

process convertReferenceFromMsfToDashesFasta {
    tag "${referenceFile.baseName}"
    
    output:
    path "${referenceFile.baseName}_ref_dashes.fasta", emit: referenceFasta

    script:
    """
    seqret -auto -stdout -sequence ${referenceFile} -sprotein1 -sformat1 unknown -osformat2 fasta -feature > ${referenceFile.baseName}_ref_dashes.fasta
    """
}

process convertReferenceFromMsfToPointsFasta {
    tag "${referenceFile.baseName}"

    output:
    path "${referenceFile.baseName}_ref_points.fasta", emit: referenceFasta

    script:
    """
    seqlim cnvt ${referenceFile} -infmt msf -outfmt fasta -o ${referenceFile.baseName}_ref_points.fasta
    """
}

workflow preProcessing {
    main:
        convertReferenceFromMsfToDashesFasta()
        convertReferenceFromMsfToPointsFasta()
    emit:
        referenceDashesFasta = convertReferenceFromMsfToDashesFasta.out.referenceFasta
        referencePointsFasta = convertReferenceFromMsfToPointsFasta.out.referenceFasta
}

workflow predictionWorkflow {
    take:
        tupleIterationTargetFile
        referenceFasta
    main:
        predictMuscle(tupleIterationTargetFile)
        sortPredictedAlignment(referenceFasta, predictMuscle.out.alignment)
        predictClustalOmega(tupleIterationTargetFile)
        predictMafft(tupleIterationTargetFile)
        predictTCoffee(tupleIterationTargetFile)
        predictHmmer(tupleIterationTargetFile, referenceFasta)
        predictKAlign(tupleIterationTargetFile)

    emit:
        muscleAlignment  = sortPredictedAlignment.out.alignment
        muscleTiming     = predictMuscle.out.timing

        clustalAlignment = predictClustalOmega.out.alignment
        clustalTiming    = predictClustalOmega.out.timing

        mafftAlignment   = predictMafft.out.alignment
        mafftTiming      = predictMafft.out.timing

        tCoffeeAlignment = predictTCoffee.out.alignment
        tCoffeeTiming    = predictTCoffee.out.timing

        hmmerAlignment   = predictHmmer.out.alignment
        hmmerTiming      = predictHmmer.out.timing

        kalignAlignment  = predictKAlign.out.alignment
        kalignTiming     = predictKAlign.out.timing
}

workflow analysisWorkflow {
    take:
        referenceFasta
        muscleAlignment
        clustalAlignment
        mafftAlignment
        tCoffeeAlignment
        hmmerAlignment
        kalignAlignment

    main:
        computeBaliScoreMuscle(muscleAlignment)
        computeBaliScoreClustal(clustalAlignment)
        computeBaliScoreMafft(mafftAlignment)
        computeBaliScoreTCoffee(tCoffeeAlignment)
        computeBaliScoreHmmer(hmmerAlignment)
        computeBaliScoreKalign(kalignAlignment)

        computeMumsaOverlapScore(
            referenceFasta,
            muscleAlignment,
            clustalAlignment,
            mafftAlignment,
            tCoffeeAlignment,
            hmmerAlignment,
            kalignAlignment
        )

        exportAlignments(
            muscleAlignment,
            clustalAlignment,
            mafftAlignment,
            tCoffeeAlignment,
            hmmerAlignment,
            kalignAlignment
        )

    emit:
        muscleBaliscore  = computeBaliScoreMuscle.out.score
        clustalBaliscore = computeBaliScoreClustal.out.score
        mafftBaliscore   = computeBaliScoreMafft.out.score
        tCoffeeBaliscore = computeBaliScoreTCoffee.out.score
        hmmerBaliscore   = computeBaliScoreHmmer.out.score
        kalignBaliscore  = computeBaliScoreKalign.out.score
}

workflow {
    preProcessing()

    tupleIterationTargetFile = Channel.from( 1 .. params.iterations )
        .map { it -> [ iteration: it, targetFile: targetFile ] }
    
    predictionWorkflow(tupleIterationTargetFile, preProcessing.out.referencePointsFasta)

    muscleAlignmentSample  = predictionWorkflow.out.muscleAlignment.randomSample(1, 1234)
    clustalAlignmentSample = predictionWorkflow.out.clustalAlignment.randomSample(1, 2345)
    mafftAlignmentSample   = predictionWorkflow.out.mafftAlignment.randomSample(1, 3456)
    tCoffeeAlignmentSample = predictionWorkflow.out.tCoffeeAlignment.randomSample(1, 4567)
    hmmerAlignmentSample   = predictionWorkflow.out.hmmerAlignment.randomSample(1, 5678)
    kalignAlignmentSample  = predictionWorkflow.out.kalignAlignment.randomSample(1, 6789)

    analysisWorkflow(
        preProcessing.out.referenceDashesFasta,
        muscleAlignmentSample,
        clustalAlignmentSample,
        mafftAlignmentSample,
        tCoffeeAlignmentSample,
        hmmerAlignmentSample,
        kalignAlignmentSample
    )
}

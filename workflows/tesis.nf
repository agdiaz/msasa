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
    predictMsasaSingleMs;
    predictMsasaBlosum;
    sortPredictedAlignment;
    predictMsasaSingleMatching;
} from "./modules/predictors/main.nf"

include {
    computeBaliScore as computeBaliScoreMuscle;
    computeBaliScore as computeBaliScoreClustal;
    computeBaliScore as computeBaliScoreMafft;
    computeBaliScore as computeBaliScoreTCoffee;
    computeBaliScore as computeBaliScoreHmmer;
    computeBaliScore as computeBaliScoreKalign;
    computeBaliScore as computeBaliScoreMsasaSingleMS;
    computeBaliScore as computeBaliScoreMsasaBlosum;
    computeBaliScore as computeBaliScoreMsasaSingleMatching;
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
        predictMsasaSingleMs(tupleIterationTargetFile)
        predictMsasaBlosum(tupleIterationTargetFile)
        predictMsasaSingleMatching(tupleIterationTargetFile)

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

        msasaSingleMsAlignment = predictMsasaSingleMs.out.alignment
        msasaSingleMsTiming    = predictMsasaSingleMs.out.timing

        msasaBlosumAlignment = predictMsasaBlosum.out.alignment
        msasaBlosumTiming    = predictMsasaBlosum.out.timing

        msasaSingleMatchingAlignment = predictMsasaSingleMatching.out.alignment
        msasaSingleMatchingTiming    = predictMsasaSingleMatching.out.timing
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
        msasaSingleMsAlignment
        msasaBlosumAlignment
        msasaSingleMatchingAlignment

    main:
        computeBaliScoreMuscle("muscle", muscleAlignment)
        computeBaliScoreClustal("clustalo", clustalAlignment)
        computeBaliScoreMafft("mafft", mafftAlignment)
        computeBaliScoreTCoffee("tcoffee", tCoffeeAlignment)
        computeBaliScoreHmmer("hmmer", hmmerAlignment)
        computeBaliScoreKalign("kalign", kalignAlignment)
        computeBaliScoreMsasaSingleMS("msasa_singlems", msasaSingleMsAlignment)
        computeBaliScoreMsasaBlosum("msasa_blosum", msasaBlosumAlignment)
        computeBaliScoreMsasaSingleMatching("msasa_singlematching", msasaSingleMatchingAlignment)

        computeMumsaOverlapScore(
            referenceFasta,
            muscleAlignment,
            clustalAlignment,
            mafftAlignment,
            tCoffeeAlignment,
            hmmerAlignment,
            kalignAlignment,
            msasaSingleMsAlignment,
            msasaBlosumAlignment,
            msasaSingleMatchingAlignment
        )

        exportAlignments(
            muscleAlignment,
            clustalAlignment,
            mafftAlignment,
            tCoffeeAlignment,
            hmmerAlignment,
            kalignAlignment,
            msasaSingleMsAlignment,
            msasaBlosumAlignment,
            msasaSingleMatchingAlignment
        )
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
    msasaSingleMsAlignmentSample  = predictionWorkflow.out.msasaSingleMsAlignment.randomSample(1, 7890)
    msasaBlosumAlignmentSample  = predictionWorkflow.out.msasaBlosumAlignment.randomSample(1, 8901)
    msasaSingleMatchingAlignmentSample  = predictionWorkflow.out.msasaSingleMatchingAlignment.randomSample(1, 9012)

    analysisWorkflow(
        preProcessing.out.referenceDashesFasta,
        muscleAlignmentSample,
        clustalAlignmentSample,
        mafftAlignmentSample,
        tCoffeeAlignmentSample,
        hmmerAlignmentSample,
        kalignAlignmentSample,
        msasaSingleMsAlignmentSample,
        msasaBlosumAlignmentSample,
        msasaSingleMatchingAlignmentSample
    )
}

#! /usr/bin/env python
import sys
import os
import argparse
import subprocess
import datetime
from Bio import SeqIO
import pysam
import operator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


#def runTCRpipe(fasta, bed, output, bam, unmapped, mapping, bases, strand, reconstruction, aaF , numIterations, thresholdScore, minOverlap,
                 #rsem, bowtie2, singleCell, path, subpath, sumF, lowQ, singleEnd, fastq, trimmomatic, transInd):
def runTCRpipe(genome, output, bam, unmapped, bases, strand, numIterations,thresholdScore, minOverlap, rsem, bowtie2, singleCell, path, sumF, lowQ, samtools, top, byExp, readOverlap, oneSide):
    checkParameters(genome, strand, singleCell, path, sumF)
    if singleCell == True:
        # TODO: Fix this, won't work for SE
        #runSingleCell(fasta, bed, output, bam, unmapped, mapping, bases, strand, reconstruction, aaF , numIterations, thresholdScore, minOverlap,
        #          rsem, bowtie2, lowQ, singleEnd, fastq, trimmomatic, transInd)
        sys.exit(0)
    if path == './':
        path = os.getcwd()
    if not path.endswith('/'):
        path = path + '/'
    finalStatDict = dict()
    tcrFout = open(sumF + '.TCRs.txt','w')
    opened = False
    for cellFolder in os.listdir(path):
        fullPath = path + cellFolder + '/'
        if((os.path.exists(fullPath)) & (os.path.isdir(fullPath))):
            sys.stdout.write(str(datetime.datetime.now()) + " Working on: " + cellFolder + '\n')
            sys.stdout.flush()
            (found, nbam, nunmapped, noutput) = formatFiles(fullPath, bam, unmapped, output)
            if not found:
                sys.stderr.write(str(datetime.datetime.now()) + " There is not a bam or unmapped file in "
                                                                "this folder, moving to the next folder\n")
                sys.stderr.flush()
            else:
                currFolder = os.path.abspath(os.path.dirname(sys.argv[0])) + '/'
                reconstruction = currFolder + '/vdj.alignment'
                if genome == 'hg38':
                    fasta = currFolder + 'Data/hg38/hg38.TCR.fa'
                    bed = currFolder + 'Data/hg38/hg38.TCR.bed'
                    mapping = currFolder + 'Data/hg38/hg38.id.name.mapping.TCR.txt'
                    aaF = currFolder + 'Data/hg38/hg38.TCR.conserved.AA.txt'
                if genome == 'mm10':
                    fasta = currFolder + 'Data/mm10/mm10.TCR.fa'
                    bed = currFolder + 'Data/mm10/mm10.TCR.bed'
                    mapping = currFolder + 'Data/mm10/mm10.gene.id.mapping.TCR.txt'
                    aaF = currFolder + 'Data/mm10/mm10.conserved.AA.txt'
                if genome == 'mm10_ncbi':
                    fasta = currFolder + 'Data/mm10_ncbi/mm10.TCR.fa'
                    bed = currFolder + 'Data/mm10_ncbi/mm10.TCR.bed'
                    mapping = currFolder + 'Data/mm10_ncbi/mm10.gene.id.mapping.TCR.txt'
                    aaF = currFolder + 'Data/mm10_ncbi/mm10.conserved.AA.txt'
                if genome == 'hg19':
                    fasta = currFolder + 'Data/hg19/hg19.TCR.fa'
                    bed = currFolder + 'Data/hg19/hg19.TCR.bed'
                    mapping = currFolder + 'Data/hg19/hg19.gene.id.mapping.TCR.txt'
                    aaF = currFolder + 'Data/hg19/hg19.conserved.AA.txt'

                runSingleCell(fasta, bed, noutput, nbam, nunmapped, mapping, bases, strand, reconstruction, aaF , numIterations, thresholdScore,
                            minOverlap, rsem, bowtie2, lowQ, samtools, top, byExp, readOverlap, oneSide)
                opened = addCellToTCRsum(cellFolder, noutput, opened, tcrFout)
                finalStatDict = addToStatDict(noutput, cellFolder, finalStatDict)
    sumFout = open(sumF + '.summary.txt','w')
    sumFout.write('sample\talpha\tbeta\n')
    for cell in sorted(finalStatDict):
        fout = cell + '\t' + finalStatDict[cell]['alpha'] + '\t' + finalStatDict[cell]['beta'] + '\n'
        sumFout.write(fout)
    sumFout.close()


def addCellToTCRsum(cellFolder, noutput, opened, tcrFout):
    if os.path.isfile(noutput + '.summary.txt'):
        currOut = open(noutput + '.summary.txt','r')
        if not opened:
            opened = True
            head = currOut.readline()
            head = 'cell\t' + head
            tcrFout.write(head)
        else:
            currOut.readline()
        l = currOut.readline()
        while l != '':
            newL = cellFolder + '\t' + l
            tcrFout.write(newL)
            l = currOut.readline()
        currOut.close()
    return opened

def addToStatDict(noutput, cellFolder, finalStatDict):
    if cellFolder in finalStatDict:
        print "Error! %s appear more than once in final stat dictionary" % cellFolder
    finalStatDict[cellFolder] = {'alpha':'Failed - found V and J segments but wasn\'t able to extend them',
                                 'beta':'Failed - found V and J segments but wasn\'t able to extend them'}
    if os.path.isfile(noutput + '.summary.txt'):
        currOut = open(noutput + '.summary.txt','r')
        msgA = 'None'
        msgB = 'None'
        currOut.readline()
        l = currOut.readline()
        while l != '':
            lArr = l.strip('\n').split('\t')
            chain = lArr[0]
            stat = lArr[1]
            if stat == 'Productive':
                if chain == 'alpha':
                    msgA = 'Productive'
                else:
                    msgB = 'Productive'
            elif stat == 'Productive (no 118 PHE found)':
                if chain == 'alpha':
                    msgA = 'Productive (no 118 PHE found)'
                else:
                    msgB = 'Productive (no 118 PHE found)'
            elif stat.startswith('Unproductive'):
                if chain == 'alpha':
                    if msgA != 'Productive':
                        msgA = 'Unproductive'
                else:
                    if msgB != 'Productive':
                        msgB = 'Unproductive'
            elif stat.startswith('Failed reconstruction'):
                if stat == 'Failed reconstruction - reached maximum number of iterations':
                    if chain == 'alpha':
                        if msgA == 'None':
                            msgA = 'Failed - reconstruction didn\'t converge'
                    else:
                        if msgB == 'None':
                            msgB = 'Failed - reconstruction didn\'t converge'
                elif stat == 'Failed reconstruction - V and J segment do not overlap':
                    if chain == 'alpha':
                        if msgA == 'None':
                            msgA = 'Failed - V and J reconstruction don\'t overlap'
                    else:
                        if msgB == 'None':
                            msgB = 'Failed - V and J reconstruction don\'t overlap'
            l = currOut.readline()
        currOut.close()
        if msgA == 'None':
            alphaJunc = noutput + '.alpha.junctions.txt'
            if (os.path.isfile(alphaJunc) == True):
                if os.stat(alphaJunc).st_size == 0:
                    msgA = 'Failed - didn\'t find any V and J segments in original mapping'
                else:
                    msgA = 'Failed - found V and J segments but wasn\'t able to extend them'
            else:
                msgA = 'Failed - didn\'t find any V and J segments in original mapping'
        if msgB == 'None':
            betaJunc = noutput + '.beta.junctions.txt'
            if (os.path.isfile(betaJunc) == True):
                if os.stat(betaJunc).st_size == 0:
                    msgB = 'Failed - didn\'t find any V and J segments in original mapping'
                else:
                    msgB = 'Failed - found V and J segments but wasn\'t able to extend them'
            else:
                msgB = 'Failed - didn\'t find any V and J segments in original mapping'

    else:
        betaJunc = noutput + '.beta.junctions.txt'
        alphaJunc = noutput + '.alpha.junctions.txt'
        if os.path.isfile(betaJunc) == True:
            if os.stat(betaJunc).st_size == 0:
                msgB = 'Failed - didn\'t find any V and J segments in original mapping'
        else:
            msgB = 'Failed - didn\'t find any V and J segments in original mapping'
        if (os.path.isfile(alphaJunc) == True):
            if os.stat(alphaJunc).st_size == 0:
                msgA = 'Failed - didn\'t find any V and J segments in original mapping'
        else:
            msgA = 'Failed - didn\'t find any V and J segments in original mapping'
    finalStatDict[cellFolder]['alpha'] = msgA
    finalStatDict[cellFolder]['beta'] = msgB
    return finalStatDict


def formatFiles(fullPath, bam, unmapped, output):
    found = True
    nbam = fullPath + bam
    if bam.startswith('/'):
        nbam = fullPath + bam[1:]
    elif bam.startswith('./'):
        nbam = fullPath + bam[2:]
    nunmapped = fullPath + unmapped
    if unmapped.startswith('/'):
        nunmapped = fullPath + unmapped[1:]
    if unmapped.startswith('./'):
        nunmapped = fullPath + unmapped[2:]
    if ((os.path.isfile(nunmapped)) & (os.path.isfile(nbam))):
        noutput = makeOutputDir(output, fullPath)
    else:
        noutput = output
        found = False
    return (found, nbam, nunmapped, noutput)

def makeOutputDir(output, fullPath):
    noutput = output
    if output.startswith('/'):
        noutput = output[1:]
    if output.startswith('./'):
        noutput = output[2:]
    if output.endswith('/'):
        noutput = output[:-1]
    if output.find('/') != -1:
        outArr = noutput.split('/')
        currPath = fullPath
        for i in range(0,len(outArr)-1):
            currPath = currPath + outArr[i] + '/'
            if not os.path.exists(currPath):
                os.makedirs(currPath)
    noutput = fullPath + noutput
    return noutput


def runSingleCell(fasta, bed, output, bam, unmapped, mapping, bases, strand, reconstruction, aaF , numIterations, thresholdScore, minOverlap,
                  rsem, bowtie2, lowQ, samtools, top, byExp, readOverlap, oneSide):
    idNameDict = makeIdNameDict(mapping)
    fastaDict = makeFastaDict(fasta)
    vdjDict = makeVDJBedDict(bed, idNameDict)
    sys.stdout.write(str(datetime.datetime.now()) + " Pre-processing alpha chain\n")
    sys.stdout.flush()
    unDictAlpha = analyzeChain(fastaDict, vdjDict, output, bam, unmapped, idNameDict, bases, 'A', strand, lowQ, top, byExp, readOverlap)
    sys.stdout.write(str(datetime.datetime.now()) + " Pre-processing beta chain\n")
    sys.stdout.flush()
    unDictBeta = analyzeChain(fastaDict, vdjDict, output, bam, unmapped, idNameDict, bases, 'B', strand, lowQ, top, byExp, readOverlap)
    sys.stdout.write(str(datetime.datetime.now()) + " Reconstructing beta chains\n")
    sys.stdout.flush()
    subprocess.call([reconstruction, output + '.beta.mapped.and.unmapped.fa', output + '.beta.junctions.txt', output + '.reconstructed.junctions.beta.fa', str(numIterations), str(thresholdScore), str(minOverlap)])
    sys.stdout.write(str(datetime.datetime.now()) + " Reconstructing alpha chains\n")
    sys.stdout.flush()
    subprocess.call([reconstruction, output + '.alpha.mapped.and.unmapped.fa', output + '.alpha.junctions.txt',
        output + '.reconstructed.junctions.alpha.fa', str(numIterations), str(thresholdScore), str(minOverlap)])
    sys.stdout.write(str(datetime.datetime.now()) + " Creating full TCR sequencing\n")
    sys.stdout.flush()
    fullTcrFileAlpha = output + '.alpha.full.TCRs.fa'
    tcrF = output + '.reconstructed.junctions.alpha.fa'
    (cSeq, cName, cId) = getCInfo(vdjDict['Alpha']['C'][0],idNameDict, fastaDict)
    createTCRFullOutput(fastaDict, tcrF, fullTcrFileAlpha, bases, idNameDict, cSeq, cName, cId, oneSide)
    fullTcrFileBeta = output + '.beta.full.TCRs.fa'
    tcrF = output + '.reconstructed.junctions.beta.fa'
    (cSeq, cName, cId) = getCInfo(vdjDict['Beta']['C'][0],idNameDict, fastaDict)
    createTCRFullOutput(fastaDict, tcrF, fullTcrFileBeta , bases, idNameDict, cSeq, cName, cId, oneSide)
    sys.stdout.write(str(datetime.datetime.now()) + " Running RSEM to quantify expression of all possible isoforms\n")
    sys.stdout.flush()
    outDirInd = output.rfind('/')
    if outDirInd != -1:
        outDir = output[:outDirInd+1]
    else:
        outDir = os.getcwd()
    runRsem(outDir, rsem, bowtie2, fullTcrFileAlpha, fullTcrFileBeta, output, samtools)
    pickFinalIsoforms(fullTcrFileAlpha, fullTcrFileBeta, output)
    bestAlpha = output + '.alpha.full.TCRs.bestIso.fa'
    bestBeta = output + '.beta.full.TCRs.bestIso.fa'
    sys.stdout.write(str(datetime.datetime.now()) + " Finding productive CDR3\n")
    sys.stdout.flush()
    aaDict = makeAADict(aaF)
    if os.path.isfile(bestAlpha):
        fDictAlpha = findCDR3(bestAlpha, aaDict, fastaDict )
    else:
        fDictAlpha = dict()
    if os.path.isfile(bestBeta):
        fDictBeta = findCDR3(bestBeta, aaDict, fastaDict )
    else:
        fDictBeta = dict()
    betaRsemOut = output + '.beta.rsem.out.genes.results'
    alphaRsemOut = output + '.alpha.rsem.out.genes.results'
    alphaBam = output + '.alpha.rsem.out.transcript.sorted.bam'
    betaBam = output + '.beta.rsem.out.transcript.sorted.bam'
    sys.stdout.write(str(datetime.datetime.now()) + " Writing results to summary file\n")
    sys.stdout.flush()
    makeSingleCellOutputFile(fDictAlpha, fDictBeta, output, betaRsemOut, alphaRsemOut, alphaBam, betaBam, fastaDict,
                             unDictAlpha, unDictBeta, idNameDict)



def analyzeChainSingleEnd(fastq, trimmomatic, transInd, bowtie2, idNameDict, output, fastaDict, bases):
    mappedReadsDictAlpha = dict()
    mappedReadsDictBeta = dict()
    lenArr = [999,50,25]
    for currLen in lenArr:
        if currLen == 999:
            steps = ['none']
        else:
            steps = ['left','right']
        for side in steps:
            if side == 'left':
                crop = 'CROP:' + str(currLen)
            elif side == 'right':
                crop = 'HEADCROP:' + str(currLen)
            else:
                crop = ''
        # TODO: make sure we delete those files
            trimFq = fastq + '.' + str(currLen) + '.' + str(side) + '.trimmed.fq'
            # TODO: use bowtie trimmer instead
            # TODO: make sure about minus strand alignment
            if crop == '':
                subprocess.call(['java','-jar', trimmomatic, 'SE','-phred33',fastq ,trimFq, 'LEADING:15','TRAILING:15', 'MINLEN:20'])
            else:
                subprocess.call(['java','-jar', trimmomatic, 'SE','-phred33',fastq ,trimFq, 'LEADING:15','TRAILING:15', crop, 'MINLEN:20'])
            samF = trimFq + '.sam'
            if bowtie2 != '':
                if bowtie2.endswith('/'):
                    bowtieCall = bowtie2 + 'bowtie2'
                else:
                    bowtieCall = bowtie2 + '/bowtie2'
            else:
                bowtieCall = 'bowtie2'
            subprocess.call([bowtieCall ,'-q --phred33  --score-min L,0,0', '-x',transInd,'-U',trimFq,'-S',samF])
            if os.path.isfile(samF):
                mappedReadsDictAlpha = findReadsAndSegments(samF, mappedReadsDictAlpha, idNameDict,'A')
                mappedReadsDictBeta = findReadsAndSegments(samF, mappedReadsDictBeta, idNameDict,'B')
    alphaOut = output + '.alpha.junctions.txt'
    alphaOutReads = output + '.alpha.mapped.and.unmapped.fa'
    betaOutReads = output + '.beta.mapped.and.unmapped.fa'
    betaOut = output + '.beta.junctions.txt'
    writeJunctionFileSE(mappedReadsDictAlpha, idNameDict, alphaOut, fastaDict, bases, 'alpha')
    writeJunctionFileSE(mappedReadsDictBeta, idNameDict, betaOut, fastaDict, bases, 'beta')
    writeReadsFileSE(mappedReadsDictAlpha, alphaOutReads, fastq)
    writeReadsFileSE(mappedReadsDictBeta, betaOutReads, fastq)
    sys.exit(1)


def writeReadsFileSE(mappedReadsDict, outReads, fastq):
    if fastq.endswith('.gz'):
        subprocess.call(['gunzip', fastq])
        newFq = fastq.replace('.gz','')
    else:
        newFq = fastq
    out = open(outReads, 'w')
    fqF = open(newFq, 'rU')
    for record in SeqIO.parse(fqF, 'fastq'):
        if record.id in mappedReadsDict:
            newRec = SeqRecord(record.seq, id = record.id, description = '')
            SeqIO.write(newRec,out,'fasta')
    out.close()
    fqF.close()
    if fastq.endswith('.gz'):
        subprocess.call(['gzip',newFq])



def writeJunctionFileSE(mappedReadsDict,idNameDict, output, fastaDict, bases, chain):
    out = open(output, 'w')
    vSegs = []
    jSegs = []
    cSegs = []
    for read in mappedReadsDict:
        for seg in mappedReadsDict[read]:
            if idNameDict[seg].find('V') != -1:
                if seg not in vSegs:
                    vSegs.append(seg)
            elif idNameDict[seg].find('J') != -1:
                if seg not in jSegs:
                    jSegs.append(seg)
            elif idNameDict[seg].find('C') != -1:
                if seg not in cSegs:
                    cSegs.append(seg)
            else:
                print "Error! not V/J/C in fasta dict"
    if len(vSegs) == 0:
        print "Did not find any V segments for " + chain + " chain"
    else:
        if len(cSegs) == 0:
            print "Did not find any C segments for " + chain + " chain"
            cSegs = ['NA']
        if len(jSegs) == 0:
            print "Did not find any J segments for " + chain + " chain"
            jSegs = ['NA']
        for vSeg in vSegs:
            for jSeg in jSegs:
                for cSeg in cSegs:
                    addSegmentToJunctionFileSE(vSeg,jSeg,cSeg,out,fastaDict, bases, idNameDict)
    out.close()






def addSegmentToJunctionFileSE(vSeg,jSeg,cSeg,out,fastaDict, bases, idNameDict):
    vSeq = fastaDict[vSeg]
    if jSeg != 'NA':
        jName = idNameDict[jSeg]
        jSeq = fastaDict[jSeg]
    else:
        jSeq = ''
        jName = 'NoJ'
    if cSeg != 'NA':
        cName = idNameDict[cSeg]
        cSeq = fastaDict[cSeg]
    else:
        cName = 'NoC'
        cSeq = ''
    jcSeq = jSeq + cSeq
    lenSeg = min(len(vSeq),len(jcSeq))
    if bases != -10:
        if lenSeg < bases:
            sys.stdout.write(str(datetime.datetime.now()) + ' Bases parameter is bigger than the length of the V or J segment, taking the length' \
                    'of the V/J segment instead, which is: ' + str(lenSeg) + '\n')
            sys.stdout.flush()
        else:
            lenSeg = bases
    jTrim = jcSeq[:lenSeg]
    vTrim = vSeq[-1*lenSeg:]
    junc = vTrim + jTrim
    recordName = vSeg + '.' + jSeg + '.' + cSeg + '(' + idNameDict[vSeg] + '-' + jName + '-' + cName + ')'
    record = SeqRecord(Seq(junc,IUPAC.ambiguous_dna), id = recordName, description = '')
    SeqIO.write(record,out,'fasta')


def findReadsAndSegments(samF, mappedReadsDict, idNameDict,chain):
    samFile = pysam.AlignmentFile(samF,'r')
    readsIter = samFile.fetch(until_eof = True)
    for read in readsIter:
        if read.is_unmapped == False:
            seg = samFile.getrname(read.reference_id)
            if seg in idNameDict:
                if idNameDict[seg].find(chain) != -1:
                    readName = read.query_name
                    if readName not in mappedReadsDict:
                        mappedReadsDict[readName] = []
                    if seg not in mappedReadsDict[readName]:
                        mappedReadsDict[readName].append(seg)
    samFile.close()
    return mappedReadsDict










def makeSingleCellOutputFile(alphaDict, betaDict, output, betaRsem, alphaRsem, alphaBam, betaBam, fastaDict,
                             unDictAlpha, unDictBeta, idNameDict):
    outF = open(output + '.summary.txt', 'w')
    outF.write('Chain\tStatus\tRank of TCR\tV\tJ\tC\tCDR3 NT\tCDR3 AA\t#reads in TCR\t#reads in CDR3\t#reads in V\t#reads in J\t#reads in C\t%unmapped reads used in the reconstruction\t# unmapped reads used in the reconstruction\t%unmapped reads in CDR3\t#unmapped reads in CDR3\tV ID\tJ ID\tC ID\n')
    if (len(alphaDict) > 0):
        writeChain(outF, 'alpha',alphaDict,alphaRsem, alphaBam, fastaDict,unDictAlpha, output, idNameDict)
    if (len(betaDict) > 0):
        writeChain(outF,'beta',betaDict, betaRsem, betaBam, fastaDict, unDictBeta, output, idNameDict)
    outF.close()


def writeChain(outF, chain,cdrDict,rsemF, bamF, fastaDict, unDict, output, idNameDict):
    writtenArr = []
    if os.path.exists(rsemF):
        noRsem = False
        (rsemDict,unRsemDict) = makeRsemDict(rsemF, cdrDict)
    else:
        noRsem = True
    for tcr in cdrDict:
        jStart = -1
        cdrInd = -1
        cInd = -1
        if cdrDict[tcr]['stat'] == 'Productive':
            isProd = True
        else:
            isProd = False
        fLine = chain + '\t' + cdrDict[tcr]['stat'] + '\t'
        if noRsem:
            rank = 'NA'
        else:
            rank = getRank(tcr, rsemDict, unRsemDict, isProd, noRsem)
        fLine += str(rank) + '\t'
        nameArr = tcr.split('.')
        fLine += nameArr[0] + '\t' + nameArr[1] + '\t' + nameArr[2] + '\t'
        fLine += cdrDict[tcr]['CDR3 NT'] + '\t' + cdrDict[tcr]['CDR3 AA'] + '\t'
        fullSeq = cdrDict[tcr]['Full Seq'].upper()
        if not noRsem:
            totalCount = findCountsInRegion(bamF, 0, len(fullSeq), tcr)
            fLine += str(totalCount) + '\t'
            cName = nameArr[5]
            while cName.endswith('_2'):
                cName = cName[:-2]
            cSeq = fastaDict[cName].upper()
            cInd = fullSeq.find(cSeq)
            if cInd == -1:
                sys.stderr.write(str(datetime.datetime.now()) + 'Error! could not find C segment sequence in the full sequence\n')
                sys.stderr.flush()
                cCounts = 'NA'
            else:
                cCounts = findCountsInRegion(bamF, cInd, len(fullSeq), tcr)
            if cdrDict[tcr]['CDR3 NT'] != 'NA':
                cdrInd = fullSeq.find(cdrDict[tcr]['CDR3 NT'].upper())
            else:
                cdrInd = -1
            if ((cdrInd == -1) & (cdrDict[tcr]['CDR3 NT'] != 'NA')):
                sys.stderr.write(str(datetime.datetime.now()) + ' Error! Cound not find CDR3 NT sequence in the full sequence\n')
                sys.stderr.flush()
            if cdrInd != -1:
                cdrCounts = findCountsInRegion(bamF, cdrInd, cdrInd + len(cdrDict[tcr]['CDR3 NT']), tcr)
                jStart = cdrInd + len(cdrDict[tcr]['CDR3 NT'])
                if cInd != -1:
                    jCounts = findCountsInRegion(bamF, jStart, cInd, tcr)
                else:
                    jCounts = findCountsInRegion(bamF, jStart, jStart+50, tcr)
                vCounts = findCountsInRegion(bamF, 0, cdrInd, tcr)
                fLine += str(cdrCounts) + '\t' + str(vCounts) + '\t' + str(jCounts) + '\t' + str(cCounts) + '\t'
            else:
                fLine += 'NA\tNA\tNA\t' + str(cCounts) + '\t'
            vId = nameArr[3]
            jId = nameArr[4]
            cId = nameArr[5]
            if cdrDict[tcr]['CDR3 NT'] != 'NA':
                (unDictRatioCDR, unCDRcount) = getUnDictRatio(bamF, cdrInd , cdrInd + len(cdrDict[tcr]['CDR3 NT']), tcr, unDict)
                (unDictRatioALL, unAllcount) = getUnDictRatio(bamF, 0 , len(fullSeq), tcr, unDict)
                fLine += str(unDictRatioALL) + '\t' + str(unAllcount) + '\t' + str(unDictRatioCDR) + '\t' + str(unCDRcount) + '\t'
            else:
                fLine += 'NA\tNA\tNA\tNA\t'
            writtenArr.append(vId)
            writtenArr.append(jId)
            fLine += vId + '\t' + jId + '\t' + cId + '\n'
        else:
            fLine += 'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t' + nameArr[3] + '\t' + nameArr[4] + '\t' + nameArr[5] + '\n'
            #print fLine

        outF.write(str(fLine))
    writeFailedReconstructions(outF, chain, writtenArr, output, idNameDict, fastaDict )

def writeFailedReconstructions(outF, chain, writtenArr, output, idNameDict, fastaDict):
    recF = output + '.reconstructed.junctions.' + chain + '.fa'
    if os.path.isfile(recF):
        f = open(recF, 'rU')
        segDict = dict()
        for tcrRecord in SeqIO.parse(f, 'fasta'):
            tcrSeq = str(tcrRecord.seq)
            if tcrSeq.find('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN') != -1:
                status = 'Failed reconstruction - reached maximum number of iterations'
                segDict = addSegmentsToDict(segDict, status, writtenArr, tcrRecord, idNameDict, fastaDict)
            elif tcrSeq.find('NNNN') != -1:
                status = 'Failed reconstruction - V and J segment do not overlap'
                segDict = addSegmentsToDict(segDict, status, writtenArr, tcrRecord, idNameDict, fastaDict)
        f.close()
        if len(segDict) > 0:
            writeSegDict(segDict, outF, chain)


def writeSegDict(segDict, outF, chain):
    for seg in segDict:
        currDict = segDict[seg]
        pairs = ''
        for pair in currDict['pairs']:
            pairs += pair + '.'
        pairs = pairs[:-1]
        if currDict['len'] > 0:
            fLine = chain + '\t' + currDict['status'] + '\t'
            rank = findCurrRank(segDict, seg, currDict['len'])
            fLine += str(rank) + '\t'
            if currDict['type'] == 'V':
                fLine += currDict['name'] + '\t' + 'paired with: ' + pairs + '\t'
            else:
                fLine +=  'paired with: ' + pairs + '\t' + currDict['name'] + '\t'
            fLine += 'NA\t' + currDict['seq'] + '\tNA\tNA\tNA\t'
            fLine += 'NA\tNA\tNA\tNA\tNA\tNA\tNA\t'
            if currDict['type'] == 'V':
                fLine += seg + '\tNA\tNA\n'
            else:
                fLine += 'NA\t' + seg + '\tNA\n'
        #print fLine
            outF.write(str(fLine))


def findCurrRank(segDict, seg, currLen):
    rank = 1
    for s in segDict:
        if s != seg:
            if segDict[s]['len'] > currLen:
                rank += 1
    return rank

def addSegmentsToDict(segDict, status, writtenArr, tcrRecord, idNameDict, fastaDict):
    head = tcrRecord.id
    headArr = head.split('.')
    vId = headArr[0]
    jId = headArr[1].split('(')[0]
    currSeqArr = tcrRecord.seq.split('N')
    vSeq = currSeqArr[0]
    jSeq = currSeqArr[-1]
    minLen = min(len(fastaDict[vId]),len(fastaDict[jId]))
    tupArr = [(vId, vSeq),(jId, jSeq)]
    for i in range(0,len(tupArr)):
        (id,seq) = tupArr[i]
        if id not in writtenArr:
            if id in segDict:
                if i == 0:
                    if idNameDict[jId] not in segDict[id]['pairs']:
                        segDict[id]['pairs'].append(idNameDict[jId])
                    if str(segDict[id]['seq'][-20:]) != str(seq[-20:]):
                        sys.stderr.write(str(datetime.datetime.now()) + ' Error! reconstructed two different sequences from the same V-segment %s\n' % id)
                        sys.stderr.flush()
                else:
                    if idNameDict[vId] not in segDict[id]['pairs']:
                        segDict[id]['pairs'].append(idNameDict[vId])
                    if str(segDict[id]['seq'][:20]) != str(seq[:20]):
                        sys.stderr.write(str(datetime.datetime.now()) + ' Error! reconstructed two different sequences from the same J-segment %s\n' % id)
                        sys.stderr.flush()
            else:
                segDict[id] = dict()
                segDict[id]['status'] = status
                segDict[id]['seq'] = seq
                segDict[id]['len'] = len(seq) - minLen
                segDict[id]['pairs'] = []

                if i == 0:
                    segDict[id]['type'] = 'V'
                    segDict[id]['pairs'].append(idNameDict[jId])
                else:
                    segDict[id]['type'] = 'J'
                    segDict[id]['pairs'].append(idNameDict[vId])
                segDict[id]['name'] = idNameDict[id]

    return segDict






def getRank(tcr, rsemDict, unRsemDict, isProd, noRsem):
    if isProd:
        currDict = rsemDict
    else:
        currDict = unRsemDict
    if not noRsem:
        currCount = currDict[tcr]
        rank = 1
        for rec in currDict:
            if rec != tcr:
                if unRsemDict[rec] > currCount:
                    rank += 1
        return rank
    else:
        return 'NA'


def getUnDictRatio(bamF, start, end, tcr, unDict):
    unMappedCount = 0
    usedArr = []
    mappedFile = pysam.AlignmentFile(bamF,"rb")
    readsIter = mappedFile.fetch(tcr, start, end)
    for read in readsIter:
        if read.is_read1 :
            newName =  read.query_name + '_1'
        else:
            newName = read.query_name + '_2'
        if newName not in usedArr:
            usedArr.append(newName)
            if newName in unDict:
                unMappedCount += 1
    mappedFile.close()
    return (float(float(unMappedCount)/len(unDict)), unMappedCount)


def findCountsInRegion(bamF, start, end, tcr):
    readsArr = []
    mappedFile = pysam.AlignmentFile(bamF,"rb")
    readsIter = mappedFile.fetch(tcr, start, end)
    for read in readsIter:
        if read.is_read1 :
            newName =  read.query_name + '_1'
        else:
            newName = read.query_name + '_2'
        if newName not in readsArr:
            readsArr.append(newName)
    mappedFile.close()
    counts = len(readsArr)
    return counts

def makeRsemDict(rsemF, cdrDict):
    fDict = dict()
    unDict = dict()
    f = open(rsemF,'r')
    f.readline()
    l = f.readline()
    while l != '':
        lArr = l.strip('\n').split('\t')
        name = lArr[1]
        if name in cdrDict:
            if cdrDict[name]['stat'] == 'Productive':
                fDict[name] = float(lArr[4])
            unDict[name] = float(lArr[4])
        l = f.readline()
    f.close()
    return (fDict,unDict)



def findCDR3(fasta, aaDict, vdjFaDict):
    f = open(fasta, 'rU')
    fDict = dict()
    for record in SeqIO.parse(f, 'fasta'):
        if record.id in fDict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! same name for two fasta entries %s\n' % record.id)
            sys.stderr.flush()
        else:
            idArr = record.id.split('.')
            vSeg = idArr[0]
            jSeg = idArr[1]
            if ((vSeg in aaDict) & (jSeg in aaDict)):
                currDict = findVandJaaMap(aaDict[vSeg],aaDict[jSeg],record.seq)
            else:
                if vSeg in aaDict:
                    newVseg = aaDict[vSeg]
                else:
                    vId = idArr[3]
                    currSeq = vdjFaDict[vId]
                    newVseg = getBestVaa(Seq(currSeq))
                if jSeg in aaDict:
                    newJseg = aaDict[jSeg]
                else:
                    jId = idArr[4]
                    currSeq= vdjFaDict[jId]
                    newJseg = getBestJaa(Seq(currSeq))
                currDict = findVandJaaMap(newVseg,newJseg,record.seq)
            fDict[record.id] = currDict
    f.close()
    return fDict



def getBestJaa(currSeq):
    firstSeq = getNTseq(currSeq)
    secondSeq = getNTseq(currSeq[1:])
    thirdSeq = getNTseq(currSeq[2:])
    pos = 10
    seq = ''
    found = False
    for s in [firstSeq, secondSeq, thirdSeq]:
        tempSeq = s[:8]
        indF = tempSeq.find('F')
        if indF != -1:
            if ((indF+3) <= len(s)):
                if s[indF+3] == 'G':
                    found = True
                    if indF < pos:
                        pos = indF
                        seq = s[indF:]
        indG = tempSeq.find('G')
        if (indG != -1):
            if found == False:
                if ((indG + 2) <= len(s)):
                    if s[indG+2] == 'G':
                        if indG < pos:
                            found = True
                            seq = s[indG:]
                            pos = indG
    if ((found == False) & (indF != -1)):
        seq = s[indF:]
    if seq != '':
        return seq
    else:
        return firstSeq

def getBestVaa(currSeq):
    firstSeq = getNTseq(currSeq)
    secondSeq = getNTseq(currSeq[1:])
    thirdSeq = getNTseq(currSeq[2:])
    pos = 10
    seq = ''
    for s in [firstSeq, secondSeq, thirdSeq]:
        #print "S: " + s
        tempSeq = s[-8:]
        #print "tempSeq: " + tempSeq
        ind = tempSeq.find('C')
        stopInd = tempSeq.find('*')
        #print "Ind: " + str(ind)
        if ((ind != -1) & (stopInd == -1)):
            #print "inside the ind"
            if ind < pos:
                goodRF = isGoodRF(s)
                if goodRF:
                    pos = ind
                    seq = s[:-8+ind + 1]
    if seq != '':
        return seq
    else:
        return firstSeq

def isGoodRF(s):
    mInd = s.find('M')
    if mInd == -1:
        return False
    stopInd = s.find('*')
    if stopInd == -1:
        return True
    stopIndNext = s[stopInd+1:].find('*')
    while stopIndNext != -1:
        stopInd = stopIndNext + stopInd + 1
        stopIndNext = s[stopInd+1:].find('*')
        mInd = s[stopInd+1:].find('M')
        mInd = mInd + stopInd + 1
    if mInd != -1:
        return True
    else:
        return False



def findVandJaaMap(vSeg,jSeg,fullSeq):
    fDict = dict()
    firstSeq = getNTseq(fullSeq)
    secondSeq = getNTseq(fullSeq[1:])
    thirdSeq = getNTseq(fullSeq[2:])
    ntArr = [fullSeq, fullSeq[1:],fullSeq[2:]]
    aaSeqsArr = [firstSeq, secondSeq, thirdSeq]
    cdrArr = []
    posArr = []
    fPharr = []
    for aaSeq in aaSeqsArr:
        (cdr, pos, curPh) = getCDR3(aaSeq, vSeg,jSeg)
        cdrArr.append(cdr)
        posArr.append(pos)
        fPharr.append(curPh)
    maxLen = 0
    bestCDR = ''
    bestSeq = ''
    hasStop = False
    bestPos = -1
    bestCDRnt = ''
    foundGood = False
    vPos = -1
    jPos = -1
    fPh = False
    for i in range(0,3):
        if posArr[i] != -1:
            if ((cdrArr[i] != 'Only J') & (cdrArr[i] != 'Only V')):
                if len(cdrArr[i]) > maxLen:
                    if cdrArr[i].find('*') == -1:
                        foundGood = True
                        bestCDR = cdrArr[i]
                        bestPos = posArr[i]
                        maxLen = len(cdrArr[i])
                        bestSeq = ntArr[i]
                        fPh = fPharr[i]
                    else:
                        if maxLen == 0:
                            foundGood = True
                            bestPos = posArr[i]
                            bestCDR = cdrArr[i]
                            maxLen = len(cdrArr[i])
                            bestSeq = ntArr[i]
                            hasStop = True
                            fPh = fPharr[i]
                else:
                    if hasStop == True:
                        if cdrArr[i].find('*') == -1:
                            foundGood = True
                            bestPos = posArr[i]
                            hasStop = False
                            bestCDR = cdrArr[i]
                            maxLen = len(cdrArr[i])
                            bestSeq = ntArr[i]
                            fPh = fPharr[i]
            else:
                if not foundGood:
                    fPh = fPharr[i]
                    if (cdrArr[i] == 'Only J'):
                        jPos = posArr[i]-i
                    elif (cdrArr[i] == 'Only V'):
                        vPos = posArr[i]-i
    if ((vPos != -1) & (jPos != -1) & (not foundGood)):
        bestCDRnt = fullSeq[3*vPos:3*jPos]
        bestCDR = 'NA'
    elif bestPos != -1:
        bestCDRnt = bestSeq[3*bestPos : 3*bestPos+3*len(bestCDR)]
    if bestCDR.find('*') != -1:
        stat = 'Unproductive - stop codon'
    elif fPh:
        stat = 'Productive (no 118 PHE found)'
    else:
        stat = 'Productive'
    if maxLen == 0:
        if (('Only J' in cdrArr) & ('Only V' in cdrArr)):
            stat = 'Unproductive - Frame shift'
        else:
            if (('Only J' not in cdrArr) & ('Only V' in cdrArr)):
                stat = 'Unproductive - found only V segment'
            elif (('Only J' in cdrArr) & ('Only V' not in cdrArr)):
                stat = 'Unproductive - found only J segment'
            elif (('Only J' not in cdrArr) & ('Only V' not in cdrArr)):
                stat = 'Unproductive - didn\'t find V and J segment'
            else:
                stat = 'Unproductive'
            bestCDR = 'NA'
            bestCDRnt = 'NA'
    fDict['stat'] = stat
    fDict['CDR3 AA'] = bestCDR
    fDict['CDR3 NT'] = bestCDRnt
    fDict['Full Seq'] = fullSeq
    return fDict



def getCDR3(aaSeq, vSeq, jSeq):
    minDist = 14
    pos = -1
    for i in range(0,len(aaSeq) - len(vSeq) + 1):
        subAA = aaSeq[i:i+len(vSeq)]
        if len(subAA) != len(vSeq):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Wrong sub length\n')
            sys.stderr.flush()
        dist = 0
        for k in range(0,len(vSeq)):
            if vSeq[k] != subAA[k]:
                dist += 1
        if ((dist < minDist) & (subAA.endswith('C'))):
            minDist = dist
            pos = i + len(vSeq)
    jPos = -1
    minDistJ = 4
    curPh = False
    for j in range(pos+1, len(aaSeq) - len(jSeq) + 1):
        subAA = aaSeq[j: j + len(jSeq)]
        if len(subAA) != len(jSeq):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Wrong subj length\n')
            sys.stderr.flush()
        dist = 0
        for m in range(0,len(jSeq)):
            if jSeq[m] != subAA[m]:
                dist += 1
        if (dist <= minDistJ):
            if isLegal(subAA):
                jPos = j
                minDistJ = dist
                curPh = False
            else:
                if dist < minDistJ:
                    curPh = True
                    jPos = j
                    minDistJ = dist
    if pos == -1:
        if jPos != -1:
            return('Only J', jPos, curPh)
        else:
            return('No V/J found', -1, curPh)
    else:
        if jPos == -1:
            return('Only V',pos, curPh)
    return(aaSeq[pos:jPos], pos, curPh )

# Checks that the conserved amino acids remain
def isLegal(subAA):
    if (len(subAA)<4):
        return False
    if (subAA[0] == 'F'):
        if ((subAA[1] == 'G') | (subAA[3] == 'G')):
            return True
    if ((subAA[1] == 'G') & (subAA[3] == 'G')):
        return True
    if ((subAA[0] == 'G') & (subAA[2] == 'G')):
        return True
    if subAA.startswith('GR'):
        return True
    if subAA.startswith('SR'):
        return True
    return False




def getNTseq(fullSeq):
    mod = len(fullSeq) % 3
    if mod != 0:
        fSeq = fullSeq[:-mod].translate()
    else:
        fSeq = fullSeq.translate()
    return fSeq

def findSeqAndLengthOfAA(aaSeq):
    fLen = 0
    fSeq = ''
    startArr = []
    stopArr = []
    startM = aaSeq.find('M')
    while startM != -1:
        startArr.append(startM)
        startM = aaSeq.find('M', startM + 1)
    stopPos = aaSeq.find('*')
    while stopPos != -1:
        stopArr.append(stopPos)
        stopPos = aaSeq.find('*', stopPos + 1)

    if ((len(startArr) == 0) | (len(stopArr) == 0)):
        return(fSeq,fLen)
    for stP in startArr:
        currStop = findStop(stP, stopArr)
        if currStop == -1:
            return (fSeq, fLen)
        else:
            currLen = currStop - stP
            if currLen >= fLen:
                fLen = currLen
                fSeq = aaSeq[stP:currStop]
    return fSeq


def findStop(stP, stopArr):
    for x in stopArr:
        if x > stP:
            return x
    return -1


def makeAADict(aaF):
    fDict = dict()
    f = open(aaF,'r')
    l = f.readline()
    while l != '':
        lArr = l.strip('\n').split('\t')
        if lArr[0] in fDict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Warning! %s appear twice in AA file\n' % lArr[0])
            sys.stderr.flush()
        fDict[lArr[0]] = lArr[1]
        l = f.readline()
    f.close()
    return fDict



def pickFinalIsoforms(fullTcrFileAlpha, fullTcrFileBeta, output):
    pickFinalIsoformChain(fullTcrFileAlpha, output + '.alpha.full.TCRs.bestIso.fa', output + '.alpha.rsem.out.genes.results')
    pickFinalIsoformChain(fullTcrFileBeta, output + '.beta.full.TCRs.bestIso.fa', output + '.beta.rsem.out.genes.results')




def pickFinalIsoformChain(fullTCRfa, newFasta, rsemF):
    if os.path.isfile(fullTCRfa):
        f = open(fullTCRfa, 'rU')
        outFa = open(newFasta, 'w')
        fastaDict = dict()
        byVJDict = dict()
        for record in SeqIO.parse(f,'fasta'):
            if record.id in fastaDict:
            #record.id = record.id + '_2'
                sys.stderr.write(str(datetime.datetime.now()) + 'error! same name for two fasta entries %s\n' % record.id)
                sys.stderr.flush()
            fastaDict[record.id] = record.seq
            onlyVJrec = str(record.id)
            idArr = onlyVJrec.strip('\n').split('.')
            vjStr = idArr[0] + '.' + idArr[1]
            if vjStr not in byVJDict:
                byVJDict[vjStr] = []
            byVJDict[vjStr].append(record.id)
        for vjStr in byVJDict:

            if len(byVJDict[vjStr]) == 1:
                cId = byVJDict[vjStr][0]
                cSeq = fastaDict[cId]
                newRec = SeqRecord(cSeq, id = cId, description = '')
                SeqIO.write(newRec,outFa,'fasta')
            else:
            #print vjStr
            #print byVJDict[vjStr]
                bestId = findBestC(byVJDict[vjStr], rsemF)
            #print "best: " + bestId
                bSeq = fastaDict[bestId]
                newRec = SeqRecord(bSeq, id = bestId, description = '')
                SeqIO.write(newRec,outFa,'fasta')
        outFa.close()
        f.close()

def findBestC(vjArr, rsemF):
    if (os.path.exists(rsemF)):
        f = open(rsemF, 'r')
        f.readline()
        l = f.readline()
        bestSeq = 'name'
        maxCount = 0.0
        while l != '':
            lArr = l.strip('\n').split('\t')
            if lArr[0] in vjArr:
                currCount = float(lArr[4])
                if currCount > maxCount:
                    bestSeq = lArr[0]
                    maxCount = currCount
            l = f.readline()
        f.close()
        if bestSeq == 'name':
            return vjArr[0]
        return bestSeq
    else:
        return vjArr[0]


def runRsem(outDir, rsem, bowtie2, fullTcrFileAlpha, fullTcrFileBeta, output, samtools):
    if samtools != '':
        if samtools[-1] != '/':
            samtools += '/'
    rsemIndDir = outDir + 'rsem_ind'
    if os.path.exists(rsemIndDir) == False:
        os.makedirs(rsemIndDir)
    if rsem != '':
        if rsem[-1] != '/':
            rsem += '/'
    if bowtie2 != '':
        if bowtie2[-1] != '/':
            bowtie2 += '/'
    if os.path.exists(fullTcrFileAlpha):
        if bowtie2 != '':
            subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2', '--bowtie2-path', bowtie2 ,
                             '-q', fullTcrFileAlpha, rsemIndDir + '/VDJ.alpha.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                             '--bowtie2', '--bowtie2-path',bowtie2, '--bowtie2-mismatch-rate', '0.0' , '--paired-end', output + '.alpha.R1.fa',
                             output + '.alpha.R2.fa', rsemIndDir + '/VDJ.alpha.seq', output + '.alpha.rsem.out'])
        else:
            subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2',
                             '-q', fullTcrFileAlpha, rsemIndDir + '/VDJ.alpha.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                             '--bowtie2', '--bowtie2-mismatch-rate', '0.0', '--paired-end', output + '.alpha.R1.fa',
                             output + '.alpha.R2.fa', rsemIndDir + '/VDJ.alpha.seq', output + '.alpha.rsem.out'])
        unsortedBam = output + '.alpha.rsem.out.transcript.bam'
        if not os.path.exists(unsortedBam):
            print "RSEM did not produce any transcript alignment files for alpha chain, please check the -rsem parameter"
        else:
            sortedBam = output + '.alpha.rsem.out.transcript.sorted.bam'
            if not os.path.exists(sortedBam):
                subprocess.call([samtools + 'samtools', 'sort','-o',sortedBam, unsortedBam])
                subprocess.call([samtools + 'samtools', 'index', sortedBam])

    else:
        sys.stdout.write(str(datetime.datetime.now()) + " Did not reconstruct any alpha chains, not running RSEM on alpha\n")
        sys.stdout.flush()
    if os.path.exists(fullTcrFileBeta):
        if bowtie2 != '':
            subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2', '--bowtie2-path', bowtie2 ,
                             '-q', fullTcrFileBeta, rsemIndDir + '/VDJ.beta.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q', '--bowtie2', '--bowtie2-path',
                             bowtie2, '--bowtie2-mismatch-rate', '0.0', '--paired-end', output + '.beta.R1.fa', output + '.beta.R2.fa',
                             rsemIndDir + '/VDJ.beta.seq', output + '.beta.rsem.out'])
        else:
            subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2',
                             '-q', fullTcrFileBeta, rsemIndDir + '/VDJ.beta.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q', '--bowtie2',
                              '--bowtie2-mismatch-rate', '0.0', '--paired-end', output + '.beta.R1.fa', output + '.beta.R2.fa',
                             rsemIndDir + '/VDJ.beta.seq', output + '.beta.rsem.out'])
        unsortedBam = output + '.beta.rsem.out.transcript.bam'
        if not os.path.exists(unsortedBam):
            print "RSEM did not produce any transcript alignment files for beta chain, please check the -rsem parameter"
        else:
            sortedBam = output + '.beta.rsem.out.transcript.sorted.bam'
            if not os.path.exists(sortedBam):
                subprocess.call([samtools + 'samtools', 'sort','-o',sortedBam, unsortedBam])
                subprocess.call([samtools + 'samtools', 'index', sortedBam])
    else:
        sys.stdout.write(str(datetime.datetime.now()) + " Did not reconstruct any beta chains, not running RSEM on beta\n")
        sys.stdout.flush()


def createTCRFullOutput(fastaDict, tcr, outName, bases, mapDict, cSeq, cName, cId, oneSide):
    tcrF = open(tcr, 'rU')
    found = False
    ffound = False
    recNameArr = []
    for tcrRecord in SeqIO.parse(tcrF, 'fasta'):
        addedC = False
        tcrSeq = str(tcrRecord.seq)
        if tcrSeq.find('NNNNN') == -1 :
            if ffound == False:
                ffound = True
                outF = open(outName, 'w')
            idArr = tcrRecord.id.split('.')
            vEns = idArr[0]
            jEns = idArr[1].split('(')[0]
            vSeq = fastaDict[vEns]
            jSeq = fastaDict[jEns]
            recNameArr = writeRecord(tcrRecord, tcrSeq, addedC, vEns, jEns, vSeq, jSeq, mapDict,bases, cSeq, cId, cName, outF,fastaDict, recNameArr)
        elif oneSide:
            curSeq = tcrSeq.split('NNNN')[0]
            jSeg = findBestJforSeq(curSeq,fastaDict,mapDict)
            if jSeg != 'NA':
                if ffound == False:
                    ffound = True
                    outF = open(outName, 'w')
                idArr = tcrRecord.id.split('.')
                vEns = idArr[0]
                vSeq = fastaDict[vEns]
                for jEns in jSeg:
                    jSeq = fastaDict[jEns]
                    recNameArr = writeRecord(tcrRecord, curSeq, addedC, vEns, jEns, vSeq, jSeq, mapDict,bases, cSeq, cId, cName, outF,fastaDict, recNameArr)
    tcrF.close()
    if found == True:
        outF.close()

def writeRecord(tcrRecord, tcrSeq, addedC, vEns, jEns, vSeq, jSeq, mapDict, bases, cSeq, cId, cName, outF,fastaDict, recNameArr):
    vSeqTrim = ''
    jSeqTrim = ''
    if bases == -10:
        bases = min(len(vSeq), len(jSeq))
    elif bases > len(jSeq):
        jSeq = jSeq + cSeq
        addedC = True
    found = False
    for i in reversed(range(20,bases)):
        juncStart = tcrSeq[:i]
        vInd = vSeq.find(juncStart)
        if (vInd != -1):
            found = True
            vSeqTrim = vSeq[:vInd]
            break
    if found == False:
        vSeqTrim = vSeq[:-bases]
    found = False
    for j in reversed(range(20,bases)):
        juncEnd = tcrSeq[-j:]
        jInd = jSeq.find(juncEnd)
        if (jInd != -1):
            found = True
            jSeqTrim = jSeq[jInd + j:]
            break
    if found == False:
        jSeqTrim = jSeq[bases:]
    # Add TRBC or TRAC
    cArr = []
    if (str(tcrRecord.id).find('TRB')!= -1):
        for ens in mapDict:
            if mapDict[ens].find('TRBC') != -1:
                cArr.append(ens)
    elif (str(tcrRecord.id).find('TRA')!= -1):
        for ens in mapDict:
            if mapDict[ens].find('TRAC') != -1:
                cArr.append(ens)
    else:
        sys.stderr.write(str(datetime.datetime.now()) + " Error! no TRBC or TRAC\n")
        sys.stderr.flush()
    if not addedC:
        for ens in cArr:
            cSeq = fastaDict[ens]
            newSeq = vSeqTrim + tcrSeq + jSeqTrim + cSeq
            newId = mapDict[vEns] + '.' + mapDict[jEns] + '.' + mapDict[ens] + '.' + vEns + '.' + jEns + '.' + ens
            while newId in recNameArr:
                newId += '_2'
            recNameArr.append(newId)
            record = SeqRecord(Seq(newSeq,IUPAC.ambiguous_dna), id = newId, description = '')
            SeqIO.write(record,outF,'fasta')
    else:
        newSeq = vSeqTrim + tcrSeq + jSeqTrim
        newId = mapDict[vEns] + '.' + mapDict[jEns] + '.' + cName + '.' + vEns + '.' + jEns + '.' + cId
        while newId in recNameArr:
            newId += '_2'
        recNameArr.append(newId)
        record = SeqRecord(Seq(newSeq,IUPAC.ambiguous_dna), id = newId, description = '')
        SeqIO.write(record,outF,'fasta')
    return recNameArr



def findBestJforSeq(curSeq,fastaDict,idNameDict):
    jArrOld = findJsPerLen(curSeq, fastaDict, idNameDict,20)
    if len(jArrOld) == 0:
        return 'NA'
    for x in range(21,len(curSeq)):
        newArr = findJsPerLen(curSeq, fastaDict, idNameDict,x)
        if len(newArr) == 0:
            return jArrOld
        else:
            jArrOld = newArr
    print 'Found a full J segment as the V/J junction, ignoring this reconstruction'
    return 'NA'

def findJsPerLen(curSeq, fastaDict, idNameDict,trim):
    fArr = []
    for seq in fastaDict:
        if idNameDict[seq].find('J') != -1:
            jSeq = fastaDict[seq]
            lenJ = len(jSeq)
            for i in range(0,lenJ):
                if ((i + trim) <= lenJ):
                    trimJ = jSeq[i:i+trim]
                    if curSeq.find(trimJ) != -1:
                        if seq not in fArr:
                            fArr.append(seq)
                            break
    return fArr




def analyzeChain(fastaDict, vdjDict, output, bam, unmapped, idNameDict, bases, chain, strand, lowQ, top, byExp, readOverlap):
    junctionSegs = makeJunctionFile(bam, chain, output, bases, vdjDict, fastaDict, idNameDict, top, byExp, readOverlap)
    unDict = writeReadsFile(bam, unmapped, junctionSegs, output, vdjDict, chain, strand, lowQ)
    return unDict

def getCInfo(bedEntry, idNameDict, fastaDict):
    bedArr = bedEntry.strip('\n').split('\t')
    cId = bedArr[3]
    cName = idNameDict[cId]
    cSeq = fastaDict[cId]
    return(cSeq,cName,cId)


def makeJunctionFile(bam, chain, output, bases, vdjDict, fastaDict, idNameDict, top, byExp, readOverlap):
    mappedFile = pysam.AlignmentFile(bam,"rb")
    if chain == 'A':
        vdjChainDict = vdjDict['Alpha']
        outName = output + '.alpha.junctions.txt'
    elif chain == 'B':
        vdjChainDict = vdjDict['Beta']
        outName = output + '.beta.junctions.txt'
    else:
        sys.stderr.write(str(datetime.datetime.now()) + ' Error! chain parameter for function analyzeChain can only be A or B\n')
        sys.stderr.flush()
    jSegs = vdjChainDict['J']
    vSegs = vdjChainDict['V']
    (cSeq, cId, cName) = getCInfo(vdjChainDict['C'][0], idNameDict, fastaDict)
    vjSegs = []
    for x in jSegs:
        vjSegs.append(x)
    for y in vSegs:
        vjSegs.append(y)
    vjReads = dict()
    (vjReads, vjCounts) = loadReadsToDict(vjSegs, mappedFile, vjReads, readOverlap)
    junctionSegs = writeJunctions(vjReads,outName, bases, fastaDict, idNameDict, cSeq, top, vjCounts, byExp)
    if len(junctionSegs) == 0:
        sys.stdout.write(str(datetime.datetime.now()) + ' Did not find any V-J reads, searching for V-C and J-C reads:\n')
        sys.stdout.flush()
        cReads = dict()
        (cReads, cCountsDict) = loadReadsToDict(vdjChainDict['C'], mappedFile, cReads, readOverlap)
        junctionSegs = writeJunctionsWithC(vjReads,outName, bases, fastaDict, idNameDict, cReads)
    mappedFile.close()
    return junctionSegs

# Similar to "writeJunctionsWithC, only that instead of looking for V-J paired-reads, it looks for
# V-C and J-C paired-reads
# INPUT:
#       vjReads - reads dict of the V and J segments created by loadReadsToDict
#       outName - output name for junction file
#       bases - number of bases to take from V and J for the junction
#       fastaDict
#       idNameDict
#       cReads - reads dict of the C segments created by loadReadsToDict
# OUTPUT:
#        fArr - the V and J segments for which we found a junction for
def writeJunctionsWithC(vjReads,outName, bases, fastaDict, idNameDict, cReads):
    out = open(outName,'w')
    fArr = []
    vArr = []
    jArr = []
    for seg in vjReads:
        if len(vjReads[seg]) > 0 :
            for cSeg in cReads:
                if (len(cReads[cSeg]) > 0) :
                    if (len([val for val in vjReads[seg]['first'] if val in cReads[cSeg]['second']]) > 0) |\
                                (len([val for val in vjReads[seg]['second'] if val in cReads[cSeg]['first']]) > 0) :
                        if idNameDict[seg].find('J') != -1 :
                            if seg not in jArr:
                                jArr.append(seg)
                        elif idNameDict[seg].find('V') != -1 :
                            if seg not in vArr:
                                vArr.append(seg)
                        fArr.append(seg)
    for vSeg in vArr:
        for jSeg in jArr:
            vSeqFa = fastaDict[vSeg]
            jSeqFa = fastaDict[jSeg]
            lenSeg = min(len(vSeqFa),len(jSeqFa))
            if bases != -10:
                if lenSeg < bases:
                    sys.stdout.write(str(datetime.datetime.now()) + ' Bases parameter is bigger than the length of the V or J segment, taking the length' \
                                        'of the V/J segment instead, which is: ' + str(lenSeg) + '\n')
                    sys.stdout.flush()
                else:
                    lenSeg = bases
            jTrim = jSeqFa[:lenSeg]
            vTrim = vSeqFa[-1*lenSeg:]
            junc = vTrim + jTrim
            recordName = vSeg + '.' + jSeg + '(' + idNameDict[vSeg] + '-' + idNameDict[jSeg] + ')'
            record = SeqRecord(Seq(junc,IUPAC.ambiguous_dna), id = recordName, description = '')
            SeqIO.write(record,out,'fasta')
    out.close()
    return fArr


# Load all the reads from a list of segments into a dictionary.
# INPUT: segsDict: A dict, where the key is a segment, and the value is an array of bed entries of this segment
#        mappedFile: Bam file of the mapped reaeds
#       readDict: A dictionary. The keys are the segment name, the values is a dictionary 'first':[] and 'second':[]
#                 where 'first' array holds the query name of R1's that overlap this segment, and 'second' holds the
#                 query name of R2's that overlap the segment.
def loadReadsToDict(segsDict, mappedFile, readDict, readOverlap):
    countDict = dict()
    for seg in segsDict:
        lArr = seg.strip('\n').split('\t')
        segName = lArr[3]
        readDict[segName] = {'first':[],'second':[]}
        lArr = seg.strip('\n').split('\t')
        chr = lArr[0]
        start = int(lArr[1])
        end = int(lArr[2])
        readsIter = mappedFile.fetch(chr, start-1, end+1)
        readCounter = 0
        for read in readsIter:
            overlap = read.get_overlap(start-1,end+1)
            if (end-start) < readOverlap:
                readOverlap = end-start-15
            if overlap >= readOverlap:
                currName = read.query_name
                if read.is_read1:
                    if currName not in readDict[segName]['first']:
                        readDict[segName]['first'].append(currName)
                        readCounter += 1
                elif read.is_read2:
                    if currName not in readDict[segName]['second']:
                        readDict[segName]['second'].append(currName)
                        readCounter += 1
        countDict[segName] = readCounter
    return (readDict, countDict)


def writeReadsFile(bam, unmapped, junctionSegs, output, vdjDict, chain, strand, lowQ):
    if chain == 'A':
        vdjChainDict = vdjDict['Alpha']
        outReads = output + '.alpha.mapped.and.unmapped.fa'
        pairedReads1 = output + '.alpha.R1.fa'
        pairedReads2 = output + '.alpha.R2.fa'
    elif chain == 'B':
        vdjChainDict = vdjDict['Beta']
        outReads = output + '.beta.mapped.and.unmapped.fa'
        pairedReads1 = output + '.beta.R1.fa'
        pairedReads2 = output + '.beta.R2.fa'
    else:
        sys.stderr.write(str(datetime.datetime.now()) + ' Error! chain parameter for function analyzeChain can only be A or B\n')
        sys.stderr.flush()
    out = open(outReads, 'w')
    constDict = vdjChainDict['C']
    # This dict for unmapped reads has reads that should be rev.comp in the revcomp arr, otherwise in id.
    # For every read, the value is a tuple - the first value is first/second, to make sure there are no errors.
    # The second value is id/revcomp, to see what sequence should be written.
    # Note: This classification is about what should happen to the unmapped reads, not how their paired maaped
    # reads were read.
    unmappedDict = dict()
    seqDict = dict()
    alignedDict = dict()
    mappedPairsDict = dict()
    lowQDict = dict()
    for seg in constDict:
        (unmappedDict, alignedDict, seqDict, mappedPairsDict, lowQDict) = addReadsToDict(unmappedDict, seg, bam, out, False, alignedDict, seqDict, strand, 'C', mappedPairsDict, lowQDict)
    vSegs = vdjChainDict['V']
    for vSeg in vSegs:
        vSegName = vSeg.strip('\n').split('\t')[3]
        if vSegName in junctionSegs:
            (unmappedDict, alignedDict, seqDict, mappedPairsDict, lowQDict) = addReadsToDict(unmappedDict, vSeg, bam, out, True, alignedDict, seqDict, strand, 'V', mappedPairsDict, lowQDict)
    jSegs = vdjChainDict['J']
    for jSeg in jSegs:
        jSegName = jSeg.strip('\n').split('\t')[3]
        if jSegName in junctionSegs:
            (unmappedDict, alignedDict, seqDict, mappedPairsDict, lowQDict) = addReadsToDict(unmappedDict, jSeg, bam, out, True, alignedDict, seqDict, strand, 'J', mappedPairsDict, lowQDict)
    unDict = dict()
    (seqDict,unDict) = writeUnmappedReads(unmappedDict, out, unmapped, seqDict, unDict, alignedDict, lowQDict, lowQ)
    seqDict = addMappedPairsToSeqDict(seqDict, bam, out, lowQ, alignedDict)
    writeSeqDict(seqDict, pairedReads1, pairedReads2)
    out.close()
    return unDict

def addMappedPairsToSeqDict(seqDict, bam, out, lowQ, alignedDict):
    firstDict = dict()
    secondDict = dict()
    for name in seqDict:
        if seqDict[name][0] == '0':
            firstDict[name] = '1'
        if seqDict[name][1] == '1':
            secondDict[name] = '1'
        if ((seqDict[name][1] == '1') & (seqDict[name][0] == '0')):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! empty record insdie seqDict\n')
            sys.stderr.flush()
    f = pysam.AlignmentFile(bam,"rb")
    readsIter = f.fetch()
    for read in readsIter:
        name = read.query_name
        pos = -1
        if ((read.query_name in firstDict) & (read.is_read1)):
            pos = 0
            check = '0'
        elif ((read.query_name in secondDict) & (read.is_read2)):
            pos = 1
            check = '1'
        if pos != -1:
            qSeq = Seq(read.query_sequence, IUPAC.ambiguous_dna)
            if read.is_read1:
                currName = name + '\\1'
            else:
                currName = name + '\\2'
            if lowQ:
                if read.mate_is_reverse:
                    if read.is_reverse:
                        rSeq = qSeq.reverse_complement()
                    else:
                        rSeq = qSeq
                else:
                    if not read.is_reverse:
                        rSeq = qSeq.reverse_complement()
                    else:
                        rSeq = qSeq
                if currName not in alignedDict:
                    alignedDict[currName] = str(rSeq)
                    record = SeqRecord(rSeq, id = currName, description = '')
                    SeqIO.write(record,out,'fasta')
                else:
                    if str(alignedDict[currName]) != str(rSeq):
                        sys.stderr.write(str(datetime.datetime.now()) + ' Error! read %s has two different sequences in alignedDict\n' % currName)
                        sys.stderr.flush()
            if read.is_reverse:
                qSeq = qSeq.reverse_complement()
            if seqDict[name][pos] != check:
                if str(qSeq) != str(seqDict[name][pos]):
                    sys.stderr.write(str(datetime.datetime.now()) + ' Error! read %s has two different mapped sequences not in the V/J region\n' % name)
                    sys.stderr.flush()
            seqDict[name][pos] = qSeq
    f.close()
    return seqDict


def writeSeqDict(seqDict, r1, r2):
    r1f = open(r1,'w')
    r2f = open(r2,'w')
    for seq in seqDict:
        if ((seqDict[seq][0] != '0') & (seqDict[seq][1]!= '1')):
            seq1 = seq
            seq2 = seq
            rec1 = SeqRecord(seqDict[seq][0], id = seq1, description = '')
            rec2 = SeqRecord(seqDict[seq][1], id = seq2, description = '')
            SeqIO.write(rec1,r1f,'fasta')
            SeqIO.write(rec2,r2f,'fasta')
        else:
            sys.stderr.write(str(datetime.datetime.now()) + ' The read %s has only one mate found, ignoring it\n' % seq)
            sys.stderr.flush()
    r1f.close()
    r2f.close()

def writeUnmappedReads(unmappedDict, out, unmapped, seqDict, unDict, alignedDict, lowQDict, lowQ):
    f = pysam.AlignmentFile(unmapped,"rb")
    readsIter = f.fetch(until_eof = True)
    for read in readsIter:
        name = read.query_name
        if name in unmappedDict:
            cName = name
            unDictName = name
            (strand , ori) = unmappedDict[name]
            if (((strand == 'first') & (read.is_read2) ) | ((strand == 'second') & (read.is_read1))):
                sys.stderr.write(str(datetime.datetime.now()) + ' Error! unmapped read is inconsistent regarding first/second read\n')
                sys.stderr.flush()
            else:
                if strand == 'first' :
                    name += '\\1'
                    unDictName += '_1'
                else:
                    name += '\\2'
                    unDictName += '_2'
                if unDictName in unDict:
                    sys.stderr.write(str(datetime.datetime.now()) + ' Error! unmapped read %s appear twice in unmapped bam file\n' % cName)
                    sys.stderr.flush()
                unDict[unDictName] = '1'
                qSeq = Seq(read.query_sequence, IUPAC.ambiguous_dna)
                if ori == 'rev':
                    qSeq = qSeq.reverse_complement()
                if name in alignedDict:
                    if alignedDict[name] != str(qSeq):
                        sys.stderr.write(str(datetime.datetime.now()) + ' Error! unmapped read %s appear twice in alignedDict with differnet seqs\n' % name)
                        sys.stderr.flush()
                else:
                    if ((name not in lowQDict) | ((name in lowQDict) & (not lowQ))):
                        alignedDict[name] = str(qSeq)
                        record = SeqRecord(qSeq, id = name, description = '')
                        SeqIO.write(record,out,'fasta')
                if ori == 'rev':
                    qSeq = qSeq.reverse_complement()
            if cName not in seqDict:
                sys.stderr.write(str(datetime.datetime.now()) + ' Error! unmapped read is in unmappedDict but not in seqDict %s\n' % read.query_name)
                sys.stderr.flush()
            else:
                if strand == 'first':
                    seqDict[cName][0] = qSeq
                else:
                    seqDict[cName][1] = qSeq
    f.close()
    return (seqDict, unDict)

# Aligned dict - all the reads (with _1/_2) that were already written to the mapped.unmapped.fa file
def addReadsToDict(unmappedDict, segBed, bam, out, mappedRead, alignedDict, seqDict, strand, segType, mappedPairsDict, lowQDict):
    bedArr = segBed.strip('\n').split('\t')
    chr = bedArr[0]
    start = int(bedArr[1])
    end = int(bedArr[2])
    mappedFile = pysam.AlignmentFile(bam,"rb")
    readsIter = mappedFile.fetch(chr, start-1, end+1)
    for read in readsIter:
        currName = read.query_name
        if currName not in seqDict:
            seqDict[currName] = ['0','1']
        if read.is_read1:
            pairRead = 'second'
            readName = currName + '\\1'
            pairName = currName + '\\2'
            seqPos = 0
        elif read.is_read2:
            pairRead = 'first'
            readName = currName + '\\2'
            pairName = currName + '\\1'
            seqPos = 1
        else:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Read is not read1 and not read2\n')
            sys.stderr.flush()
        currSeq = Seq(read.query_sequence,IUPAC.ambiguous_dna)
        if read.is_reverse:
            pairOr = 'id'
            readStrand = 'minus'
            currSeq = currSeq.reverse_complement()
        else:
            readStrand = 'plus'
            pairOr = 'rev'
        seqDict[currName][seqPos] = currSeq
        if read.mate_is_unmapped:
            takePair = toTakePair(segType, strand, readStrand)
            if takePair == False:
                lowQDict[pairName] = '1'
            #takePair = True
            #if takePair:
            if currName in unmappedDict:
                if unmappedDict[currName] != (pairRead, pairOr):
                    sys.stderr.write(str(datetime.datetime.now()) + ' Error! Read %s has more than one unmppaed mate with differnet strand/mate\n' % currName)
                    sys.stderr.flush()
            unmappedDict[currName] = (pairRead, pairOr)
        if mappedRead == True :
            if readName in alignedDict:
                if alignedDict[readName] != read.query_sequence:
                    sys.stderr.write(str(datetime.datetime.now()) + ' Error! Read %s has two instances but different seuqences\n' % read.query_name)
                    sys.stderr.flush()
            else:
                alignedDict[readName] = read.query_sequence
                record = SeqRecord(Seq(read.query_sequence, IUPAC.ambiguous_dna), id = readName, description = '')
                SeqIO.write(record,out,'fasta')
    mappedFile.close()
    #print "Removed counter: " + str(cc)
    return(unmappedDict, alignedDict, seqDict, mappedPairsDict, lowQDict)


### For minus, add an unmapped pair if the current mate is: 1. V and Plus 2. J/C and minus
### For plus, exactly the opposite

def toTakePair(segType, strand, readStrand):
    if strand == 'none':
        return True
    if ((readStrand != 'minus') & (readStrand != 'plus')):
        sys.stderr.write(str(datetime.datetime.now()) + ' Error! Read strand should be plus or minus only\n')
        sys.stderr.flush()
        return True
    if ((segType == 'C') | (segType == 'J')):
        if strand == 'minus':
            if readStrand == 'minus':
                return True
            else:
                return False
        else:
            if readStrand == 'minus':
                return False
            else:
                return True
    else:
        if strand == 'minus':
            if readStrand == 'minus':
                return False
            else:
                return True
        else:
            if readStrand == 'minus':
                return True
            else:
                return False
    return True


def writeJunctions(vjReads,outName, bases, fastaDict, idNameDict, cSeq, top, vjCountsDict, byExp):
    out = open(outName,'w')
    fArr = []
    pairCountDict = dict()
    seqsDict = dict()
    for seg in vjReads:
        if idNameDict[seg].find('J') != -1 :
            if len(vjReads[seg]) > 0 :
                for sSeg in vjReads:
                    if ((idNameDict[sSeg].find('V') != -1) & (len(vjReads[sSeg]) > 0)) :
                        if (len([val for val in vjReads[seg]['first'] if val in vjReads[sSeg]['second']]) > 0) |\
                                (len([val for val in vjReads[seg]['second'] if val in vjReads[sSeg]['first']]) > 0) :
                            if seg not in fArr:
                                fArr.append(seg)
                            if sSeg not in fArr:
                                fArr.append(sSeg)
                            vSeq = fastaDict[sSeg]
                            jSeq = fastaDict[seg]
                            lenSeg = min(len(vSeq),len(jSeq))
                            if bases != -10:
                                if lenSeg < bases:
                                    if bases > len(vSeq):
                                        sys.stdout.write(str(datetime.datetime.now()) + ' Bases parameter is bigger than the length of the V segment, taking the length' \
                                              'of the V/J segment instead, which is: ' + str(lenSeg) + '\n')
                                        sys.stdout.flush()
                                    else:
                                        sys.stdout.write(str(datetime.datetime.now()) + ' Bases parameter is bigger than the length of the J segment, appending the C segment to the J segment\n')
                                        sys.stdout.flush()
                                        jSeq = jSeq + cSeq
                                        lenSeg = bases
                                else:
                                    lenSeg = bases
                            jTrim = jSeq[:lenSeg]
                            vTrim = vSeq[-1*lenSeg:]
                            junc = vTrim + jTrim
                            recordName = sSeg + '.' + seg + '(' + idNameDict[sSeg] + '-' + idNameDict[seg] + ')'
                            record = SeqRecord(Seq(junc,IUPAC.ambiguous_dna), id = recordName, description = '')
                            curCont = vjCountsDict[seg] + vjCountsDict[sSeg]
                            pairCountDict[record.seq] = curCont
                            seqsDict[record.seq] = record
    sorted_pairs = sorted(pairCountDict.items(), key=operator.itemgetter(1), reverse=True)
    if ((top == -1) | (top > len(sorted_pairs))):
        for rec in seqsDict:
            SeqIO.write(seqsDict[rec] ,out,'fasta')
    else:
        if not byExp:
            for i in range(0,top):
                SeqIO.write(seqsDict[sorted_pairs[i][0]],out,'fasta')
        else:
            wrote = 1
            SeqIO.write(seqsDict[sorted_pairs[0][0]],out,'fasta')
            curCount = sorted_pairs[0][1]
            wroteSecond = False
            for i in range(1,len(sorted_pairs)):
                if sorted_pairs[i][1] == curCount:
                    if not wroteSecond:
                        wroteSecond = True
                        SeqIO.write(seqsDict[sorted_pairs[i][0]],out,'fasta')
                        wrote += 1
                else:
                    curCount = sorted_pairs[i][1]
                    wroteSecond = False
                    SeqIO.write(seqsDict[sorted_pairs[i][0]],out,'fasta')
                    wrote += 1
                if wrote == top:
                    break

    out.close()
    return fArr





# Create a dict {'Alpha':{'C':[bed],'V':[bed],'J':[bed]}, 'Beta':{'C':[],'V':[],'J':[]}}
def makeVDJBedDict(bed,idNameDict):
    fDict = {'Alpha':{'C':[],'V':[],'J':[]}, 'Beta':{'C':[],'V':[],'J':[]}}
    f = open(bed, 'r')
    l = f.readline()
    while l != '':
        lArr = l.strip('\n').split('\t')
        gID = lArr[3]
        gName = idNameDict[gID]
        chain = ''
        if (gName.startswith('TRA')):
            chain = 'Alpha'
        elif (gName.startswith('TRB')):
            chain = 'Beta'
        else:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! %s name is not alpha or beta chain, ignoring it\n' % gName)
            sys.stderr.flush()
        if gName.find('C') != -1:
            fDict[chain]['C'].append(l)
        elif gName.find('V') != -1:
            fDict[chain]['V'].append(l)
        elif gName.find('J') != -1:
            fDict[chain]['J'].append(l)
        l = f.readline()
    f.close()
    return fDict





# Creates a dictionary of ENSEMBL ID -> fasta sequence
def makeFastaDict(fasta):
    inF = open(fasta,'rU')
    fastaDict = dict()
    for record in SeqIO.parse(inF, 'fasta'):
        fastaDict[record.id] = str(record.seq)
    inF.close()
    return fastaDict


# Creates a dictionary of ENSEMBL ID -> Gene name
def makeIdNameDict(mapping):
    f = open(mapping, 'r')
    fDict = dict()
    linesArr = f.read().split('\n')
    f.close()
    for line in linesArr:
        lineArr = line.split('\t')
        id = lineArr[0]
        name = lineArr[1]
        ind = name.find('Gene:')
        if ind != -1:
            name = name[ind+len('Gene:'):]
        if id in fDict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! %s appear twice in mapping file\n' % id)
            sys.stderr.flush()
        fDict[id] = name
    return fDict


def checkParameters(genome, strand, singleCell, path, sumF):
    if ((genome != 'hg38') & (genome != 'mm10') & (genome != 'hg19') & (genome != 'mm10_ncbi')):
        sys.exit("-genome only accept one of the following: mm10, mm10_ncbi, hg38, hg19")
    if strand.lower() not in ['none','minus','plus']:
        sys.exit("-strand should be one of: none, minus, plus")
    if not singleCell:
        if path == '':
            sys.exit("when running on multiple cells you must include the -path parameter")
        if sumF == '':
            sys.exit("when running on multiple cells you must include the -sumF parameter")
        if not os.path.isdir(path):
            sys.exit("%s path does not exists. Please check your -path parameter and run again" % path)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-genome','-g','-G', help='Alignment genome. Currently supported: mm10, mm10_ncbi and hg38', required=True)
    parser.add_argument('-singleCell', help='add if you are only running on a single cell. If so,'
                                                        'it will ignore -path and -subpath arguments', action='store_true')
    parser.add_argument('-lowQ', help='add if you want to add \"low quality\" reads as input to the reconstruction '
                                                        'algorithm', action='store_true')
    parser.add_argument('-oneSide', help='add if you want to observe reconstrctuion only from the V side', action='store_true')
    parser.add_argument('-path','-p','-P', help='The path for the data directory. Assumes that every subdirectory'
                                                    'is a single cell', default='')
    parser.add_argument('-sumF', help='prefix for summary outputs', default='')
    parser.add_argument('-bowtie2','-bw','-BW', help='Path to bowtie2. If not used assumes that bowtie2 is in the'
                                                 'default path', default = '')
    parser.add_argument('-rsem','-RSEM', help='Path to rsem. If not used assumes that rsem is in the'
                                                'default path', default = '/data/yosef/users/safik/bin/rsem-1.2.21/')
    parser.add_argument('-strand', help='Strand of the right most read in genomic coordinates. Options are: [minus, plus, '
                                        'none]. Defualt is minus', default = 'minus')
    parser.add_argument('-output','-out','-o','-O', help='output prefix, relative to /path/singleCellFolder', required=True)
    parser.add_argument('-bam', help='Input bam alignment file, relative to /path/singleCellFolder/ if working on multiple files', default = './tophat_output/picard_output/sorted.bam')
    parser.add_argument('-unmapped','-u','-U', help='bam file of the unmapped reads, relative to /path/singleCellFolder/', default = './tophat_output/unmapped.bam')
    parser.add_argument('-bases','-b','-B', help='Number of bases to take from each V and J segments, default is min(len(V), len(J) ', type=int, default=-10)
    parser.add_argument('-iterations','-iter','-i','-I', help='Number of iterations for the reconstruction'
                                                              'algorithm, default is 20', type=int, default=20)
    parser.add_argument('-samtools', help='Path to samtools. If not used assumes that samtools is in the default path', default = '')
    parser.add_argument('-score','-sc','-SC', help='Alignment score threshold. Default is 15', type=int, default=15)
    parser.add_argument('-top','-t','-T', help='Take only the top x combination of V and J, based on the sum '
                                               'number of reads that map to both. Default is to take all', type=int, default=-1)
    parser.add_argument('-readOverlap','-ro','-readoverlap', help='Add a read to list of mapped reads only if it maps at least X bases'
                                               'to the V/J/C segment. Default is 1', type=int, default=1)
    parser.add_argument('-byExp', help='if using the Top option, add this tag if you want to take only two chains from each'\
                                                        'read count, until top is reached', action='store_true')

    parser.add_argument('-overlap','-ol','-OL', help='Number of minimum bases that overlaps V and J ends,'
                                                              'default is 10', type=int, default=10)
    args = parser.parse_args()
    runTCRpipe(args.genome, args.output, args.bam, args.unmapped, args.bases, args.strand,
                args.iterations,args.score, args.overlap, args.rsem, args.bowtie2,
                  args.singleCell, args.path, args.sumF, args.lowQ, args.samtools, args.top, args.byExp, args.readOverlap,
                  args.oneSide)


DINARA_PHASING_DEBUG_READ=695 \
/home/kokyriakidis/Downloads/dinara-build/Executable/dinara \
--input /home/kokyriakidis/Downloads/GIAB_HG002_TESTS/PAW70337.chr1_35_45Mb_subregion.fasta \
--Reads.minReadLength 1000 \
--threads 20 \
--Kmers.k 50 \
--Assembly.mode3.minAnchorCoverage 4 \
--Assembly.mode3.maxAnchorCoverage 400 \
--Align.minAlignedMarkerCount 2 \
--Align.minAlignedFraction 0 \
--OverlapCandidates.method InvertedIndex \
--OverlapCandidates.driftRateTolerance 0.05 \
--OverlapCandidates.maxEndFuzz 0 \
--OverlapCandidates.minOverlapLength 0 \
--MarkerGraph.alwaysSave \
--memoryMode filesystem \
--memoryBacking disk \
--assemblyDirectory /home/kokyriakidis/Downloads/PAW70337.chr1_35_45Mb_subregion.fasta


/home/kokyriakidis/Downloads/dinara-build/Executable/dinara \
--input /home/kokyriakidis/Downloads/GIAB_HG002_TESTS/PAW70337.HIC2_4Mb.subregion.fasta \
--config Nanopore-r10.4.1_e8.2-400bps_sup-Herro-Jan2025 \
--Reads.minReadLength 1000 \
--threads 20 --Kmers.k 50 \
--Kmers.minimizerW 50 \
--Kmers.minMarkerSpanFraction 0.0 \
--Assembly.mode3.maxAnchorCoverage 400 \
--OverlapCandidates.method InvertedIndex \
--OverlapCandidates.driftRateTolerance 0.05 \
--OverlapCandidates.maxEndFuzz 0 \
--Align.minAlignedMarkerCount 2 \
--OverlapCandidates.minOverlapLength 0 \
--Assembly.mode3.minAnchorCoverage 2 \
--Assembly.mode3.minWindowBaseSpan 10000 \
--Assembly.mode3.minInterWindowCoverage 4 \
--Assembly.mode3.minInterWindowEdgeCoverage 1 \
--Assembly.mode3.minCommonForBackbone 0 \
--MarkerGraph.alwaysSave \
--memoryMode filesystem \
--memoryBacking disk \
--assemblyDirectory /home/kokyriakidis/Downloads/PAW70337.HIC2_4Mb.subregion.fasta 




micromamba activate samtools

HG002v1_1_REF="/home/kokyriakidis/Downloads/hg002v1.1.fasta"

samtools faidx ${HG002v1_1_REF} chr12_PATERNAL \
> chr12_PATERNAL.fasta

samtools faidx \
/home/kokyriakidis/Downloads/Assembly.fasta \
"10-25-0-0-P0" \
> /home/kokyriakidis/Downloads/chr12_PATERNAL_10-25-0-0-P0.fasta


micromamba activate minimap2                                                                   

ASSEMBLY_SEGMENT_PATERNAL="/home/kokyriakidis/Downloads/chr12_PATERNAL_10-25-0-0-P0.fasta"

minimap2 \
-t 20 \
-I 20G \
-cx asm10 \
--eqx \
"chr12_PATERNAL.fasta" \
"chr12_PATERNAL_10-25-0-0-P0.fasta" \
> /home/kokyriakidis/Downloads/chr12_PATERNAL_10-25-0-0-P0.paf


/home/kokyriakidis/Downloads/DisplayPafAlignments \
--reference chr12_PATERNAL.fasta \
--assembly chr12_PATERNAL_10-25-0-0-P0.fasta \
--paf /home/kokyriakidis/Downloads/chr12_PATERNAL_10-25-0-0-P0.paf \
--html /home/kokyriakidis/Downloads/chr12_PATERNAL_10-25-0-0-P0_diplayPafAlignments.html










/home/kokyriakidis/Downloads/dinara-build/Executable/dinara \
--input /home/kokyriakidis/Downloads/GIAB_HG002_TESTS/GIAB_HG002_PAW70337_RAW_chr1_15-15.4.fastq \
--config Nanopore-r10.4.1_e8.2-400bps_sup-Herro-Jan2025 \
--Reads.minReadLength 1000 \
--threads 20 --Kmers.k 50 \
--Kmers.minimizerW 50 \
--Kmers.minMarkerSpanFraction 0.0 \
--Assembly.mode3.maxAnchorCoverage 200 \
--OverlapCandidates.method InvertedIndex \
--OverlapCandidates.driftRateTolerance 0.05 \
--OverlapCandidates.maxEndFuzz 0 \
--Align.minAlignedMarkerCount 2 \
--OverlapCandidates.minOverlapLength 1000 \
--Assembly.mode3.minAnchorCoverage 2 \
--Assembly.mode3.minWindowBaseSpan 0 \
--Assembly.mode3.minInterWindowCoverage 3 \
--Assembly.mode3.minInterWindowEdgeCoverage 0 \
--Assembly.mode3.minCommonForBackbone 0 \
--Assembly.mode3.hetDropHomopolymer 1 \
--Assembly.mode3.hetDropRepeat 1 \
--Assembly.mode3.hetMinVaf 0.12 \
--Assembly.mode3.hetMinSupport 6 \
--MarkerGraph.alwaysSave \
--memoryMode filesystem \
--memoryBacking disk \
--assemblyDirectory /home/kokyriakidis/Downloads/GIAB_HG002_PAW70337_RAW_chr1_15-15.4.fastq




/home/kokyriakidis/Downloads/dinara-build/Executable/dinara \
--input /home/kokyriakidis/Downloads/E821_StdMix_TESTS/E821_StdMix_RAW_chr1-15-30.fasta  \
--config Nanopore-r10.4.1_e8.2-400bps_sup-Herro-Jan2025 \
--Reads.minReadLength 1000 \
--threads 20 \
--Kmers.k 50 \
--Kmers.minimizerW 50 \
--Kmers.minMarkerSpanFraction 0.0 \
--Assembly.mode3.maxAnchorCoverage 200 \
--OverlapCandidates.method InvertedIndex \
--OverlapCandidates.driftRateTolerance 0.05 \
--OverlapCandidates.maxEndFuzz 0 \
--Align.minAlignedMarkerCount 2 \
--OverlapCandidates.minOverlapLength 1000 \
--Assembly.mode3.minAnchorCoverage 2 \
--Assembly.mode3.minWindowBaseSpan 5000 \
--Assembly.mode3.minInterWindowCoverage 3 \
--Assembly.mode3.minInterWindowEdgeCoverage 0 \
--Assembly.mode3.minCommonForBackbone 0 \
--Assembly.mode3.hetDropHomopolymer 0 \
--Assembly.mode3.hetDropRepeat 0 \
--Assembly.mode3.hetMinVaf 0.20 \
--Assembly.mode3.hetMinSupport 6 \
--MarkerGraph.alwaysSave \
--memoryMode filesystem \
--memoryBacking disk \
--assemblyDirectory /home/kokyriakidis/Downloads/E821_StdMix_RAW_chr1-15-30.fasta






/home/kokyriakidis/Downloads/E821_StdMix_TESTS/E821_StdMix_RAW_chr1-15-30.fasta 

/home/kokyriakidis/Downloads/shasta2-build/Executable/shasta2 \
--config /home/kokyriakidis/Downloads/shasta2-build/shasta2-install/conf/ONT-diploid-40x-60kb-Q25-2026.06.26.conf \
--min-read-length 1000 \
--transitive-reduction-max-distance 0 \
--threads 20 \
--memory-mode filesystem \
--memory-backing disk \
--keep-binary-data \
--k 2 \
--read-following-segment-length-threshold 5000 \
--min-anchor-coverage 2 \
--max-anchor-coverage 500 \
--min-anchor-graph-edge-coverage 1 \
--write-intermediate-assembly-stages \
--external-anchors-name /home/kokyriakidis/Downloads/GIAB_HG002_PAW70337_RAW_chr1_15-15.4.fastq/Shasta2ExternalAnchors  \
--external-anchor-graph-name /home/kokyriakidis/Downloads/GIAB_HG002_PAW70337_RAW_chr1_15-15.4.fastq/Data/Shasta2-Shasta2AnchorGraph \
--input /home/kokyriakidis/Downloads/E821_StdMix_TESTS/E821_StdMix_RAW_chr1-15-30.fasta \
--output /home/kokyriakidis/Downloads/E821_StdMix_TESTS/E821_StdMix_RAW_chr1-15-30.fasta




./shasta2-build/Executable/shasta2 \
--config /home/kokyriakidis/Downloads/shasta2-build/shasta2-install/conf/ONT-diploid-40x-60kb-Q25-2026.06.26.conf \
--min-read-length 1000 \
--threads 20 \
--memory-mode filesystem \
--memory-backing disk \
--keep-binary-data \
--write-intermediate-assembly-stages \
--input /home/kokyriakidis/Downloads/E821_StdMix_TESTS/E821_StdMix_RAW_chr1-15-30.fasta \
--output /home/kokyriakidis/Downloads/shasta2_raw_E821_StdMix_RAW_chr1-15-30.fasta

/home/kokyriakidis/Downloads/shasta2-build/Executable/shasta2 \
--config /home/kokyriakidis/Downloads/shasta2-build/shasta2-install/conf/ONT-diploid-40x-60kb-Q25-2026.06.26.conf \
--min-read-length 1000 \
--transitive-reduction-max-distance 0 \
--threads 20 \
--memory-mode filesystem \
--memory-backing disk \
--keep-binary-data \
--k 2 \
--read-following-segment-length-threshold 5000 \
--min-anchor-coverage 2 \
--max-anchor-coverage 500 \
--min-anchor-graph-edge-coverage 1 \
--write-intermediate-assembly-stages \
--external-anchors-name /home/kokyriakidis/Downloads/E821_StdMix_RAW_chr1-15-30.fasta/Shasta2ExternalAnchors  \
--external-anchor-graph-name /home/kokyriakidis/Downloads/E821_StdMix_RAW_chr1-15-30.fasta/Data/Shasta2-Shasta2AnchorGraph \
--input /home/kokyriakidis/Downloads/E821_StdMix_TESTS/E821_StdMix_RAW_chr1-15-30.fasta \
--output /home/kokyriakidis/Downloads/shasta2_E821_StdMix_RAW_chr1-15-30.fasta




/home/kokyriakidis/Downloads/shasta2-build/Executable/shasta2 \
--config /home/kokyriakidis/Downloads/shasta2-build/shasta2-install/conf/ONT-diploid-40x-60kb-Q25-2026.06.26.conf \
--min-read-length 1000 \
--transitive-reduction-max-distance 0 \
--threads 20 \
--memory-mode filesystem \
--memory-backing disk \
--keep-binary-data \
--k 2 \
--read-following-segment-length-threshold 10000 \
--min-anchor-coverage 4 \
--max-anchor-coverage 500 \
--min-anchor-graph-edge-coverage 4 \
--write-intermediate-assembly-stages \
--external-anchors-name /home/kokyriakidis/Downloads/GIAB_HG002_PAW70337_RAW_chr1_15-15.4.fastq/Shasta2ExternalAnchors  \
--external-anchor-graph-name /home/kokyriakidis/Downloads/GIAB_HG002_PAW70337_RAW_chr1_15-15.4.fastq/Data/Shasta2-Shasta2AnchorGraph \
--input /home/kokyriakidis/Downloads/GIAB_HG002_TESTS/GIAB_HG002_PAW70337_RAW_chr1_15-15.4.fastq  \
--output /home/kokyriakidis/Downloads/shasta2_GIAB_HG002_PAW70337_RAW_chr1_15-15.4.fastq



/home/kokyriakidis/Downloads/shasta2-build/Executable/shasta2 \
--config /home/kokyriakidis/Downloads/shasta2-build/shasta2-install/conf/ONT-diploid-40x-60kb-Q25-2026.06.26.conf \
--min-read-length 1000 \
--threads 20 \
--memory-mode filesystem \
--memory-backing disk \
--keep-binary-data \
--write-intermediate-assembly-stages \
--input /home/kokyriakidis/Downloads/GIAB_HG002_TESTS/GIAB_HG002_PAW70337_RAW_chr1_15-15.4.fastq

./shasta2-build/Executable/shasta2 \
--config ./conf/ONT-diploid-40x-60kb-Q25-2026.04.15.conf \
--min-read-length 1000 \
--transitive-reduction-max-distance 0 \
--threads 20 \
--memory-mode filesystem \
--memory-backing disk \
--keep-binary-data \
--k 50 \
--read-following-segment-length-threshold 10000 \
--min-anchor-coverage 4 \
--max-anchor-coverage 500 \
--min-anchor-graph-edge-coverage 4 \
--write-intermediate-assembly-stages \
--external-anchors-name /home/kokyriakidis/Downloads/PAW70337.HIC2_4Mb.subregion.fasta/Shasta2ExternalAnchors  \
--external-anchor-graph-name /home/kokyriakidis/Downloads/PAW70337.HIC2_4Mb.subregion.fasta/Shasta2ExternalAnchorGraph  \
--input /home/kokyriakidis/Downloads/GIAB_HG002_TESTS/PAW70337.HIC2_4Mb.subregion.fasta  \
--output /home/kokyriakidis/Downloads/shasta2_PAW70337.HIC2_4Mb.subregion.fasta







/home/kokyriakidis/Downloads/shasta2-build/Executable/shasta2 \
--command explore \
--explore-access unrestricted \
--output /home/kokyriakidis/Downloads/shasta2_E821_StdMix_RAW_chr1-15-30.fasta
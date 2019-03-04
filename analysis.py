import os
from enum import Enum
import tempfile
from tempfile import TemporaryDirectory
from itertools import chain, groupby, filterfalse
from os.path import abspath, join, basename, dirname, splitext, realpath, commonprefix
from os import makedirs
import attr
from functools import wraps

from dmpy.distributedmake import DistributedMake, get_dm_arg_parser


@attr.s
class GenomeReference(object):
    fasta = attr.ib()
    index = attr.ib()
    ref_dict = attr.ib()


class VCFType(Enum):
    vcf_compressed = ("z",)

    def __init__(self, bcftools_cli_flag):
        self.bcftools_cli_flag = bcftools_cli_flag


@attr.s
class DMWorker(object):
    dm = attr.ib()
    global_tmpdir = attr.ib()
    tmpdir = attr.ib(init=False)

    bcftools_stats_targets = attr.ib(default=attr.Factory(set))

    def __attrs_post_init__(self):
        self.tmpdir = TemporaryDirectory(dir=self.global_tmpdir)

    def cleanup(self):
        self.tmpdir.cleanup()

    def unzip_to(self, unzipped, zipped):
        self.dm.add(unzipped, zipped,
                    "gzip -dc {} > {}".format(zipped, unzipped))

    def add_read_groups(self, bam_with_rg, fastq):
        self.dm.add(bam_with_rg, fastq.star_bam,
                    ("picard AddOrReplaceReadGroups"
                     " SO=coordinate"
                     " I={0} O={1} RGPL=ILLUMINA"
                     # Each sample got a single library prep,
                     # so lib and sample identifiers can be the same
                     " RGSM={2} RGLB={2}"
                     # Instrument:Run_number:Flowcell_id:lane
                     " RGID=$(samtools view {0} | head -1 | cut -f1-4 -d':')"
                     # The ID + library prep id
                     " RGPU=$(samtools view {0} | head -1 | cut -f1-4 -d':'):{2}"
                     ).format(fastq.star_bam, bam_with_rg, fastq.experimental_sample))
        return bam_with_rg

    def mark_duplicates(self, dedupped_bam, bam, remove_duplicates=False):
        command = (
            "picard MarkDuplicates"
            " I={} O={}"
            " CREATE_INDEX=true"
            " VALIDATION_STRINGENCY=SILENT"
            " M=output.metrics"
        ).format(bam, dedupped_bam)
        if remove_duplicates:
            command += " REMOVE_DUPLICATES=true"
        self.dm.add(dedupped_bam, bam, command)

    def merge_bams(self, merged_bam, bams, add_read_group=False):
        if add_read_group:
            tmpdir = tempfile.mkdtemp(dir=self.tmpdir.name)
            commands = [
                "cd {}".format(tmpdir),
                " && ",
                " && ".join(
                    ["ln -f -s {} {}".format(bam, "{}.bam".format(idx)) for idx, bam in
                     enumerate(bams)]
                ),
                " && ",
                "samtools merge -r {} *.bam".format(merged_bam),
            ]

        else:
            commands = [
                "samtools merge {}".format(merged_bam),
            ]
            for bam in bams:
                commands.append(" " + bam)

        commands.append(" && samtools index " + merged_bam)

        self.dm.add(merged_bam, bams, ''.join(commands))

    def split_n_trim(self, trimmed_bam, input_bam, genome_reference_object):
        self.dm.add(trimmed_bam,
                    [input_bam, genome_reference_object.index, genome_reference_object.ref_dict],
                    (
                        "GenomeAnalysisTK -Xmx10g"
                        " -T SplitNCigarReads "
                        "-R {} -I {} -o {} "
                        "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 "
                        "-U ALLOW_N_CIGAR_READS"
                    ).format(genome_reference_object.fasta, input_bam, trimmed_bam))

    def index_fasta(self, fasta):
        index = fasta + ".fai"
        self.dm.add(index, fasta, "samtools faidx {}".format(fasta))
        return index

    def index_fasta_for_broad(self, ensembl_reference_unzip):
        index = self.index_fasta(ensembl_reference_unzip)
        ref_dict = self.dict_fasta(ensembl_reference_unzip)
        return GenomeReference(fasta=ensembl_reference_unzip, index=index, ref_dict=ref_dict)

    def dict_fasta(self, unzipped_reference):
        ref_dict = splitext(unzipped_reference)[0] + '.dict'
        self.dm.add(ref_dict, unzipped_reference,
                    "picard CreateSequenceDictionary R={} O={}".format(unzipped_reference,
                                                                       ref_dict))
        return ref_dict

    def bqsr(self, output_bam, input_bam, known_sites, ref_fasta):
        recal_dir = output_bam + "_data"
        recal_table = join(recal_dir, "recal_data.table")
        post_recal_table = join(recal_dir, "post_recal_data.table")
        recal_plots = join(recal_dir, 'recal_comparison.pdf')

        self.bqsr_generate_recal_table(recal_table, input_bam, known_sites, ref_fasta)
        self.bqsr_generate_post_recal_table(post_recal_table, input_bam, recal_table, known_sites,
                                            ref_fasta)
        self.bqsr_generate_recal_plots(recal_plots, recal_table, post_recal_table, ref_fasta)
        self.bqsr_apply_recal_table(output_bam, input_bam, recal_table, known_sites, ref_fasta)

    def bqsr_generate_recal_table(self, output_table, input_bam, known_sites, ref_fasta):
        if isinstance(known_sites, str):
            known_sites = [known_sites]
        command = (
            'GenomeAnalysisTK -T BaseRecalibrator'
            ' -R {}'
            ' -I {}'
            ' -o {}'
        ).format(ref_fasta, input_bam, output_table)
        for known_sites_file in known_sites:
            command += ' -knownSites {}'.format(known_sites_file)
        self.dm.add(output_table, known_sites + [input_bam], command)

    def bqsr_generate_post_recal_table(self, output_table, input_bam, input_table, known_sites,
                                       ref_fasta):
        if isinstance(known_sites, str):
            known_sites = [known_sites]
        command = (
            'GenomeAnalysisTK -T BaseRecalibrator'
            ' -R {}'
            ' -I {}'
            ' -BQSR {}'
            ' -o {}'
        ).format(ref_fasta, input_bam, input_table, output_table)
        for known_sites_file in known_sites:
            command += ' -knownSites {}'.format(known_sites_file)
        self.dm.add(output_table, known_sites + [input_bam, input_table], command)

    def bqsr_generate_recal_plots(self, output_plots, before_table, after_table, ref_fasta):
        command = (
            'GenomeAnalysisTK -T AnalyzeCovariates'
            ' -R {}'
            ' -before {}'
            ' -after {}'
            ' -plots {}'
        ).format(ref_fasta, before_table, after_table, output_plots)
        self.dm.add(output_plots, [before_table, after_table], command)

    def bqsr_apply_recal_table(self, output_bam, input_bam, recal_table, known_sites, ref_fasta):
        if isinstance(known_sites, str):
            known_sites = [known_sites]
        command = (
            'GenomeAnalysisTK -T PrintReads'
            ' -R {}'
            ' -I {}'
            ' -o {}'
            ' -BQSR {}'
        ).format(ref_fasta, input_bam, output_bam, recal_table)
        self.dm.add(output_bam, known_sites + [recal_table, input_bam], command)

    def call_variants(self, output_vcf, input_bam, ref_fasta):
        command = ('GenomeAnalysisTK -T HaplotypeCaller'
                   ' -R {}'
                   ' -I {}'
                   ' -dontUseSoftClippedBases'
                   ' -stand_call_conf 20.0'
                   ' -o {}').format(ref_fasta, input_bam, output_vcf)
        self.dm.add(output_vcf, [input_bam, ref_fasta], command)

    def set_variant_filters(self, output_vcf, input_vcf, ref_fasta,
                            filters=""):
        command = ("GenomeAnalysisTK -T VariantFiltration"
                   " -R {}"
                   " -V {}"
                   ' -o {}').format(ref_fasta, input_vcf, output_vcf)
        command += filters
        self.dm.add(output_vcf, input_vcf, command)

    def apply_variant_filters(self, output_vcf, input_vcf, filters=None,
                              output_type=VCFType.vcf_compressed):
        if filters is None:
            filters = ['.', 'PASS']
        command = "bcftools view {} -o {} -O {}".format(input_vcf, output_vcf,
                                                        output_type.bcftools_cli_flag)
        command += ' -f ' + ','.join(filters)
        self.dm.add(output_vcf, input_vcf, command)

    def bcftools_stats(self, input_vcf, output_file=None):
        if output_file is None:
            output_file = input_vcf + '.stats'
        command = 'bcftools stats {} > {}'.format(input_vcf, output_file)
        self.dm.add(output_file, input_vcf, command)
        self.bcftools_stats_targets.add(output_file)


@attr.s(slots=True)
class FastqMessenger(object):
    fastq = attr.ib()
    fastq_sample = attr.ib()
    experimental_sample = attr.ib()
    fastqc_report = attr.ib(default=None)
    star_bam = attr.ib(default=None)
    star_bam_bai = attr.ib(default=None)
    star_bam_with_rg = attr.ib(default=None)
    read_dist_report = attr.ib(default=None)


def gather_fastq_info(rna_seq_dir):
    fastqs = []
    for sample_dir in os.listdir(rna_seq_dir):
        for fastq in os.listdir(join(rna_seq_dir, sample_dir)):
            fastq_full_path = join(rna_seq_dir, sample_dir, fastq)
            fastq_sample = fastq.rpartition(".fastq.gz")[0]
            fastqs.append(
                FastqMessenger(fastq=fastq_full_path,
                               fastq_sample=fastq_sample,
                               experimental_sample=sample_dir)
            )
    return fastqs


@attr.s(slots=True)
class SampleMessenger(object):
    fastqs = attr.ib()
    name = attr.ib()
    star_bam_merged = attr.ib(default=None)
    star_bam_with_rg_merged = attr.ib(default=None)
    star_bam_with_rg_merged_dedupped = attr.ib(default=None)
    star_bam_with_rg_merged_dedupped_snt = attr.ib(default=None)
    star_bam_with_rg_merged_dedupped_snt_recal = attr.ib(default=None)
    variant_calls = attr.ib(default=None)
    variant_calls_with_filters = attr.ib(default=None)
    variant_calls_filtered_general = attr.ib(default=None)
    variant_calls_filtered_for_geneiase = attr.ib(default=None)
    htseq_counts = attr.ib(default=None)
    featureCounts_counts = attr.ib(default=None)
    featureCounts_counts_simple = attr.ib(default=None)


def collate_fastqs_to_experimental_samples(fastqs):
    fastq_exp_sample = lambda x: x.experimental_sample
    sample_sorted_fastqs = sorted(fastqs, key=fastq_exp_sample)
    exp_samples = []
    for exp_samp, fastq_iter in groupby(sample_sorted_fastqs, key=fastq_exp_sample):
        exp_samples.append(
            SampleMessenger(
                fastqs=list(fastq_iter),
                name=exp_samp,
            )
        )
    return exp_samples


MAX_READ_LENGTH = 75
MAPPING_THREADS = 4
FLOW_CELL_ID = 'HYF5WBGX2'

GLOBAL_TMPDIR = abspath("tmp")
makedirs(GLOBAL_TMPDIR,exist_ok=True)

# This is the location where Kiran downloaded the FASTQ files to
RNA_SEQ_DIR = abspath("data/AdHa160517-40613574")
GENOME_ANNOTATION_DIR = abspath("data/QC/annotation")
ENSEMBL_ANNOTATIONS = abspath("data/annotations/ensembl/Homo_sapiens.GRCh38.89.gtf.gz")
ENSEMBL_REFERENCE = abspath(
    "data/reference/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
# REFSEQ_REFERENCE = abspath("data/reference/refseq/GCF_000001405.36_GRCh38.p10_genomic.fna.gz")
# REFSEQ_ANNOTATIONS = abspath(
#    "data/annotations/refseq/GCF_000001405.36_GRCh38.p10_genomic.gff.gz")
GO_TERMS = abspath("data/gene_sets/Human_GOALL_with_GO_iea_June_01_2017_symbol.gmt")

RESULTS_DIR = abspath("results")
FASTQC_DIR = join(RESULTS_DIR, "fastqc")
RSEQC_DIR = join(RESULTS_DIR, "rseqc")
ALIGNMENT_DIR = join(RESULTS_DIR, "alignment")
VARIANT_DIR = join(RESULTS_DIR, "variants")
HTSEQ_DIR = join(RESULTS_DIR, "htseq")
FEATURECOUNTS_DIR = join(RESULTS_DIR, "featureCounts")
ENSEMBL_ANNOT_DIR = join(RESULTS_DIR, "annotations", "ensembl")
ENSEMBL_REF_DIR = join(RESULTS_DIR, "reference", "ensembl")
GENE_SET_DIR = join(RESULTS_DIR, "gene_sets")

BROAD_BUNDLE_KNOWN_SITES = [
    'data/broad/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz',
    'data/broad/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
    'data/broad/bundle/hg38/dbsnp_146.hg38.vcf.gz',
]


def main():
    args = get_dm_arg_parser('RNA-seq analysis').parse_args()
    dm = DistributedMake(args_object=args, keep_going=True)
    dmw = DMWorker(dm, GLOBAL_TMPDIR)

    ensembl_reference_unzip = join(ENSEMBL_REF_DIR, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    dmw.unzip_to(ensembl_reference_unzip, ENSEMBL_REFERENCE)
    ensembl_reference_object = dmw.index_fasta_for_broad(ensembl_reference_unzip)

    ensembl_annotations_unzip = join(ENSEMBL_ANNOT_DIR, "Homo_sapiens.GRCh38.89.gtf")
    dm.add(ensembl_annotations_unzip, ENSEMBL_ANNOTATIONS,
           "gzip -dc {} > {}".format(ENSEMBL_ANNOTATIONS, ensembl_annotations_unzip))

    ensembl_annotations_real_genes = join(ENSEMBL_ANNOT_DIR,
                                          "Homo_sapiens.GRCh38.89.protein_coding.gtf")
    dm.add(ensembl_annotations_real_genes, ensembl_annotations_unzip,
           (' grep -e \'^#\' -e \'gene_biotype "protein_coding"\' {}'
            ' > {}').format(ensembl_annotations_unzip, ensembl_annotations_real_genes)
           )

    # create a new STAR index from GRCh37
    star_ref_flag = join(RESULTS_DIR, "reference/GRCh38.89_star/flag")
    star_ref_dir = dirname(star_ref_flag)
    dm.add(star_ref_flag, [ensembl_reference_unzip, ENSEMBL_ANNOTATIONS],
           (
               'STAR'
               ' --runThreadN 32'
               ' --runMode genomeGenerate'
               ' --genomeDir {0}'
               ' --genomeFastaFiles {1}'
               ' --sjdbGTFfile {2}'
               ' --sjdbOverhang {3}'
               ' && touch {4}'
           ).format(star_ref_dir,
                    ensembl_reference_unzip,
                    ensembl_annotations_unzip,
                    MAX_READ_LENGTH - 1,
                    star_ref_flag))

    # Let's first do some quality control
    fastqc_reports = []
    fastqs = gather_fastq_info(RNA_SEQ_DIR)
    for fastq in fastqs:
        fastq.fastqc_report = join(FASTQC_DIR,
                                   fastq.experimental_sample,
                                   basename(fastq.fastq).rpartition(".fastq.gz")[
                                       0] + "_fastqc.html")
        out_dir = dirname(fastq.fastqc_report)
        dm.add(fastq.fastqc_report, fastq.fastq,
               "fastqc --outdir {} {}".format(out_dir, fastq.fastq))

    # Let's map the reads using hisat2 - oh wait hisat2 quality sucks
    # STAR it is.
    for fastq in fastqs:
        fastq.star_bam = join(ALIGNMENT_DIR, "star", fastq.fastq_sample,
                              "Aligned.sortedByCoord.out.bam")
        out_dir = dirname(fastq.star_bam)
        dm.add(
            fastq.star_bam,
            [fastq.fastq, star_ref_flag],
            " ".join(
                [
                    "STAR --runThreadN {} --outSAMstrandField intronMotif".format(MAPPING_THREADS),
                    "--genomeDir", star_ref_dir,
                    "--readFilesIn", fastq.fastq,
                    "--outFileNamePrefix {}/".format(out_dir),
                    "--outSAMtype", "BAM", "SortedByCoordinate",
                    "--readFilesCommand", "zcat",
                ]
            )
        )
        fastq.star_bam_bai = fastq.star_bam + ".bai"
        dm.add(fastq.star_bam_bai, fastq.star_bam, "samtools index " + fastq.star_bam)

    exp_samples = collate_fastqs_to_experimental_samples(fastqs)

    perform_broad_variant_calling(dmw, fastqs, exp_samples, ensembl_reference_object,
                                  BROAD_BUNDLE_KNOWN_SITES, VARIANT_DIR)

    # QC summary
    multiqc_report = join(RESULTS_DIR, "multiqc", "multiqc_report.html")
    deps = list(chain(*[(f.fastqc_report, f.star_bam) for f in fastqs]))
    dm.add(multiqc_report, deps,
           "multiqc -d -dd 1 -f --outdir {} {} {}".format(dirname(multiqc_report), FASTQC_DIR,
                                                          ALIGNMENT_DIR))

    # combine experimental samples
    for sample in exp_samples:
        star_bams = [f.star_bam for f in sample.fastqs]
        sample.star_bam_merged = join(ALIGNMENT_DIR, "star_by_experimental_sample", sample.name,
                                      sample.name + ".bam")
        dmw.merge_bams(sample.star_bam_merged, star_bams, add_read_group=True)

        # try featureCounts - count reads on exons
        sample.featureCounts_counts = join(FEATURECOUNTS_DIR, sample.name, "exons",
                                           sample.name + ".exon.counts")
        dm.add(
            sample.featureCounts_counts,
            [sample.star_bam_merged, ensembl_annotations_real_genes],
            " featureCounts"
            " -t exon -F GTF -g gene_id -a {} -o {} {}".format(ensembl_annotations_real_genes,
                                                               sample.featureCounts_counts,
                                                               sample.star_bam_merged, )
        )

        # reformat counts to simpler format
        sample.featureCounts_counts_simple = sample.featureCounts_counts + ".simple_counts"
        dm.add(
            sample.featureCounts_counts_simple,
            sample.featureCounts_counts,
            "grep -v '^#' {} | cut -f 1,7 > {}".format(sample.featureCounts_counts,
                                                       sample.featureCounts_counts_simple)
        )

    sample_multiqc_report = join(RESULTS_DIR, "multiqc_per_experimental_sample",
                                 "multiqc_report.html")
    deps = [s.featureCounts_counts for s in exp_samples]
    dm.add(sample_multiqc_report,
           deps,
           "multiqc -f --outdir {} {}".format(dirname(sample_multiqc_report), FEATURECOUNTS_DIR))

    # calculate gene lengths
    fc_gene_lengths = join(FEATURECOUNTS_DIR, "gene_lengths.tsv")
    dm.add(fc_gene_lengths, exp_samples[0].featureCounts_counts,
           ("grep '^ENSG' {}"
            " | cut -f1,3,4"
            " | perl -ane 'print $F[0] . \"\\t\" . ($F[2] - $F[1]) . \"\\n\"' > {} ").format(
               exp_samples[0].featureCounts_counts, fc_gene_lengths))

    # collate counts into one csv
    exp_samples = sorted(exp_samples, key=lambda x: x.name)
    fc_counts = [s.featureCounts_counts_simple for s in exp_samples]
    fc_merged_counts = join(FEATURECOUNTS_DIR, "all_sample_counts.csv")
    assert len(fc_counts) >= 2
    dm.add(
        fc_merged_counts, fc_counts,
        "rm -f {}".format(fc_merged_counts) +
        " && join {} {}".format(*fc_counts[0:2]) + "".join(
            [" | join - " + f for f in fc_counts[2:]]) + " >> " + fc_merged_counts
    )

    # convert go terms to go id only
    go_terms_only_symbols = join(GENE_SET_DIR,
                                 "Human_GOALL_with_GO_iea_June_01_2017_symbol.gmt")
    dm.add(go_terms_only_symbols, GO_TERMS,
           "perl -pe 's/.*\%//' < {} > {}".format(GO_TERMS, go_terms_only_symbols))

    # run R analysis
    flag = join(HTSEQ_DIR, "analysis.R.flag")
    # dm.add(flag, htseq_merged_counts, "Rscript analysis.R")

    dm.execute()
    dmw.cleanup()


def in_directory(file, directory):
    # make both absolute
    directory = realpath(directory)
    file = realpath(file)

    # return true, if the common prefix of both is equal to directory
    # e.g. /a/b/c/d.rst and directory is /a/b, the common prefix is /a/b
    return commonprefix([file, directory]) == directory


def perform_broad_variant_calling(dmw, fastqs, experimental_samples, genome_reference_object,
                                  known_sites, variant_dir):
    """Broad pipeline for variant calling on RNA-seq data

    This analysis broadly follows GENOMEANALYSISTK best practices as defined at
    https://software.broadinstitute.org/GenomeAnalysisTK/documentation/article.php?id=3891,
    retrieved 2017-08-14
    """

    for fastq in fastqs:
        fastq.star_bam_with_rg = splitext(fastq.star_bam)[0] + '.with_rg.bam'
        dmw.add_read_groups(fastq.star_bam_with_rg, fastq)
    for sample in experimental_samples:
        sample.star_bam_with_rg_merged = join(ALIGNMENT_DIR,
                                              "star_by_experimental_sample",
                                              sample.name,
                                              sample.name + ".with_rg.bam")
        bams = [f.star_bam_with_rg for f in sample.fastqs]
        dmw.merge_bams(sample.star_bam_with_rg_merged, bams)

        sample.star_bam_with_rg_merged_dedupped = \
            splitext(sample.star_bam_with_rg_merged)[0] + '.dedup.bam'
        dmw.mark_duplicates(sample.star_bam_with_rg_merged_dedupped,
                            sample.star_bam_with_rg_merged,
                            remove_duplicates=True)

        sample.star_bam_with_rg_merged_dedupped_snt = \
            splitext(sample.star_bam_with_rg_merged_dedupped)[0] + '.snt.bam'
        dmw.split_n_trim(sample.star_bam_with_rg_merged_dedupped_snt,
                         sample.star_bam_with_rg_merged_dedupped,
                         genome_reference_object)
        sample.star_bam_with_rg_merged_dedupped_snt_recal = \
            splitext(sample.star_bam_with_rg_merged_dedupped_snt)[0] + '.recal.bam'
        dmw.bqsr(sample.star_bam_with_rg_merged_dedupped_snt_recal,
                 sample.star_bam_with_rg_merged_dedupped_snt,
                 known_sites,
                 genome_reference_object.fasta)

        sample.variant_calls = join(variant_dir,
                                    "star_by_experimental_sample",
                                    sample.name,
                                    sample.name + ".vcf.gz")
        dmw.call_variants(sample.variant_calls, sample.star_bam_with_rg_merged_dedupped_snt_recal,
                          genome_reference_object.fasta)
        dmw.bcftools_stats(sample.variant_calls)

        sample.variant_calls_with_filters = join(
            dirname(sample.variant_calls),
            'with_filters',
            sample.name + '.with_filters.vcf.gz'
        )
        dmw.set_variant_filters(sample.variant_calls_with_filters,
                                sample.variant_calls,
                                genome_reference_object.fasta,
                                filters=(' -window 35 -cluster 3'
                                         ' -filterName FS -filter "FS > 30.0"'
                                         ' -filterName QD -filter "QD < 2.0"'
                                         ' -filterName DP -filter "DP < 10"'))
        dmw.bcftools_stats(sample.variant_calls_with_filters)

        sample.variant_calls_filtered_general = join(
            dirname(sample.variant_calls),
            'filtered_for_general',
            sample.name + '.filtered_for_general.vcf.gz'
        )
        dmw.apply_variant_filters(sample.variant_calls_filtered_general,
                                  sample.variant_calls_with_filters,
                                  ['.', 'PASS', 'DP'])
        dmw.bcftools_stats(sample.variant_calls_filtered_general)

        sample.variant_calls_filtered_for_geneiase = join(
            dirname(sample.variant_calls), 'filtered_for_geniase',
            sample.name + '.filtered_for_geneiase.vcf.gz'
        )
        dmw.apply_variant_filters(sample.variant_calls_filtered_for_geneiase,
                                  sample.variant_calls_with_filters,
                                  ['.', 'PASS'])
        dmw.bcftools_stats(sample.variant_calls_filtered_for_geneiase)

    variant_multiqc_report = join(variant_dir, 'multiqc_report-variants.html')
    dmw.dm.add(variant_multiqc_report, dmw.bcftools_stats_targets,
               'multiqc -m bcftools --filename {} -f {}'.format(variant_multiqc_report,
                                                                variant_dir))


if __name__ == '__main__':
    main()

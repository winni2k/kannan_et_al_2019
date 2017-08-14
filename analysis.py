import os
from tempfile import TemporaryDirectory
from itertools import chain, groupby
from os.path import abspath, join, basename, dirname

from dmpy.distributedmake import DistributedMake, get_dm_arg_parser


class FrozenClass(object):
    __isfrozen = False

    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("Cannot set {}: {} is a frozen class".format(key, self))
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True


class Messenger(FrozenClass):
    def __init__(self, **kwargs):
        self.__dict__ = kwargs
        self._freeze()


class DMWorker(object):
    def __init__(self, dm, global_tmpdir):
        self.dm = dm
        self.tmpdir = TemporaryDirectory(dir=global_tmpdir)

    def cleanup(self):
        self.tmpdir.cleanup()

    def unzip_to(self, zipped, unzipped):
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
                     " RGID=$$(samtools view {0} | head -1 | cut -f1-4 -d':')"
                     # The ID + library prep id
                     " RGPU=$$(samtools view {0} | head -1 | cut -f1-4 -d':'):{2}"
                     ).format(fastq.star_bam, bam_with_rg, fastq.experimental_sample))
        return bam_with_rg

    def mark_duplicates(self, dedupped_bam, bam):
        self.dm.add(dedupped_bam, bam, (
            "picard MarkDuplicates"
            " I={} O={}"
            " CREATE_INDEX=true"
            " VALIDATION_STRINGENCY=SILENT"
            " M=output.metrics"
        ).format(bam, dedupped_bam))

    def merge_bams(self, merged_bam, bams, add_read_group=False):
        if add_read_group:
            tmpdir = TemporaryDirectory(dir=self.tmpdir.name)
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


def gather_fastq_info(rna_seq_dir):
    fastqs = []
    for sample_dir in os.listdir(rna_seq_dir):
        for fastq in os.listdir(join(rna_seq_dir, sample_dir)):
            fastq_full_path = join(rna_seq_dir, sample_dir, fastq)
            fastq_sample = fastq.rstrip(".fastq.gz")
            fastqs.append(
                Messenger(
                    fastq=fastq_full_path,
                    fastq_sample=fastq_sample,
                    experimental_sample=sample_dir,
                    fastqc_report=None,
                    star_bam=None,
                    star_bam_bai=None,
                    star_bam_with_rg=None,
                    read_dist_report=None,
                )
            )
    return fastqs


def collate_fastqs_to_experimental_samples(fastqs):
    fastq_exp_sample = lambda x: x.experimental_sample
    sample_sorted_fastqs = sorted(fastqs, key=fastq_exp_sample)
    exp_samples = []
    for exp_samp, fastq_iter in groupby(sample_sorted_fastqs, key=fastq_exp_sample):
        exp_samples.append(
            Messenger(
                fastqs=list(fastq_iter),
                name=exp_samp,
                star_bam_merged=None,
                star_bam_with_rg_merged=None,
                star_bam_with_rg_merged_dedupped=None,
                htseq_counts=None,
                featureCounts_counts=None,
                featureCounts_counts_simple=None,
            )
        )
    return exp_samples


MAX_READ_LENGTH = 75
MAPPING_THREADS = 4
FLOW_CELL_ID = 'HYF5WBGX2'

GLOBAL_TMPDIR = abspath("tmp")

# This is the location where Kiran downloaded the FASTQ files to
RNA_SEQ_DIR = abspath("data/AdHa160517-40613574")
GENOME_ANNOTATION_DIR = abspath("data/QC/annotation")
ENSEMBL_ANNOTATIONS = abspath("data/annotations/ensembl/Homo_sapiens.GRCh38.89.gtf.gz")
ENSEMBL_REFERENCE = abspath(
    "data/reference/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
REFSEQ_REFERENCE = abspath("data/reference/refseq/GCF_000001405.36_GRCh38.p10_genomic.fna.gz")
REFSEQ_ANNOTATIONS = abspath(
    "data/annotations/refseq/GCF_000001405.36_GRCh38.p10_genomic.gff.gz")
GO_TERMS = abspath("data/gene_sets/Human_GOALL_with_GO_iea_June_01_2017_symbol.gmt")

RESULTS_DIR = abspath("results")
FASTQC_DIR = join(RESULTS_DIR, "fastqc")
RSEQC_DIR = join(RESULTS_DIR, "rseqc")
ALIGNMENT_DIR = join(RESULTS_DIR, "alignment")
HTSEQ_DIR = join(RESULTS_DIR, "htseq")
FEATURECOUNTS_DIR = join(RESULTS_DIR, "featureCounts")
ENSEMBL_ANNOT_DIR = join(RESULTS_DIR, "annotations", "ensembl")
ENSEMBL_REF_DIR = join(RESULTS_DIR, "reference", "ensembl")
GENE_SET_DIR = join(RESULTS_DIR, "gene_sets")


def main():
    args = get_dm_arg_parser('RNA-seq analysis').parse_args()
    dm = DistributedMake(args_object=args, keep_going=True)
    dmw = DMWorker(dm, GLOBAL_TMPDIR)

    ensembl_reference_unzip = join(ENSEMBL_REF_DIR, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    dmw.unzip_to(ENSEMBL_REFERENCE, ensembl_reference_unzip)

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
                                   basename(fastq.fastq).rstrip(".fastq.gz") + "_fastqc.html")
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

    perform_variant_calling(dmw, fastqs, exp_samples)

    # QC summary
    multiqc_report = join(RESULTS_DIR, "multiqc", "multiqc_report.html")
    dm.add(multiqc_report, chain(*[(f.fastqc_report, f.star_bam) for f in fastqs]),
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
    dm.add(sample_multiqc_report,
           chain(*[(s.featureCounts_counts,) for s in exp_samples]),
           "multiqc -f --outdir {} {} {}".format(dirname(sample_multiqc_report), HTSEQ_DIR,
                                                 FEATURECOUNTS_DIR))

    # calculate gene lengths
    fc_gene_lengths = join(FEATURECOUNTS_DIR, "gene_lengths.tsv")
    dm.add(fc_gene_lengths, exp_samples[0].featureCounts_counts,
           ("grep '^ENSG' {}"
            " | cut -f1,3,4"
            " | perl -ane 'print $$F[0] . \"\\t\" . ($$F[2] - $$F[1]) . \"\\n\"' > {} ").format(
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


def perform_variant_calling(dmw, fastqs, experimental_samples):
    for fastq in fastqs:
        fastq.star_bam_with_rg = fastq.star_bam.rpartition('.bam')[0] + '.with_rg.bam'
        dmw.add_read_groups(fastq.star_bam_with_rg, fastq)
    for sample in experimental_samples:
        sample.star_bam_with_rg_merged = join(ALIGNMENT_DIR, "star_by_experimental_sample",
                                              sample.name,
                                              sample.name + ".with_rg.bam")

        bams = [f.star_bam_with_rg for f in sample.fastqs]
        dmw.merge_bams(sample.star_bam_with_rg_merged, bams)

        sample.star_bam_with_rg_merged_dedupped = \
            sample.star_bam_with_rg_merged.rpartition('.bam')[0] + '.dedup.bam'
        dmw.mark_duplicates(sample.star_bam_with_rg_merged_dedupped, sample.star_bam_with_rg_merged)


if __name__ == '__main__':
    main()

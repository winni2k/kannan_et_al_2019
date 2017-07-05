import os
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
                htseq_counts=None,
                featureCounts_counts=None,
                featureCounts_counts_simple=None,
            )
        )
    return exp_samples


args = get_dm_arg_parser('RNA-seq analysis').parse_args()
dm = DistributedMake(args_object=args, keep_going=True)

activate_str = ". activate rna_seq_py27"

# This is the location where Kiran downloaded the FASTQ files to
rna_seq_dir = abspath("data/AdHa160517-40613574")
star_genome_dir = abspath("data/reference/hg19_Gencode14.overhangt75_star")
genome_annotation_dir = abspath("data/QC/annotation")
ensembl_annotations = abspath("data/annotations/ensembl/Homo_sapiens.GRCh37.87.gtf.gz")

results_dir = abspath("results")
fastqc_dir = join(results_dir, "fastqc")
rseqc_dir = join(results_dir, "rseqc")
alignment_dir = join(results_dir, "alignment")
htseq_dir = join(results_dir, "htseq")
featureCounts_dir = join(results_dir, "featureCounts")

# create a new STAR index from GRCh37
# dm.add()

# Pre-processing of the gtf file
ensembl_annotations_with_chr = join(results_dir, "annotations", "ensembl",
                                    "Homo_sapiens.GRCh37.87.with_chr.gtf.gz")
dm.add(ensembl_annotations_with_chr, ensembl_annotations,
       (
           "gzip -dc {}"
           " | perl -pe 's/^([^#])/chr$$1/'"
           " | perl -pe 's/^chrMT/chrM/'"
           " | pigz -c > {}"
       ).format(ensembl_annotations, ensembl_annotations_with_chr))

# and provide an unzipped version
ensembl_annotations_with_chr_no_zip = join(results_dir, "annotations", "ensembl",
                                           "Homo_sapiens.GRCh37.87.with_chr.gtf")
dm.add(ensembl_annotations_with_chr_no_zip, ensembl_annotations_with_chr,
       "gzip -dc {} > {}".format(ensembl_annotations_with_chr, ensembl_annotations_with_chr_no_zip))

# Let's first do some quality control
fastqc_reports = []
fastqs = gather_fastq_info(rna_seq_dir)
for fastq in fastqs:
    fastq.fastqc_report = join(fastqc_dir,
                               fastq.experimental_sample,
                               basename(fastq.fastq).rstrip(".fastq.gz") + "_fastqc.html")
    out_dir = dirname(fastq.fastqc_report)
    dm.add(fastq.fastqc_report, fastq.fastq,
           "fastqc --outdir {} {}".format(out_dir, fastq.fastq))


def extract_experimental_sample(fastq):
    return fastq.experimental_sample


# Let's map the reads using hisat2 - oh wait hisat2 quality sucks
# STAR it is.
for fastq in fastqs:
    fastq.star_bam = join(alignment_dir, "star", fastq.fastq_sample,
                          "Aligned.sortedByCoord.out.bam")
    out_dir = dirname(fastq.star_bam)
    dm.add(
        fastq.star_bam,
        fastq.fastq,
        " ".join(
            [
                "STAR --runThreadN 2 --outSAMstrandField intronMotif",
                "--genomeDir", star_genome_dir,
                "--readFilesIn", fastq.fastq,
                "--outFileNamePrefix {}/".format(out_dir),
                "--outSAMtype", "BAM", "SortedByCoordinate",
                "--readFilesCommand", "zcat",
            ]
        )
    )
    fastq.star_bam_bai = fastq.star_bam + ".bai"
    dm.add(fastq.star_bam_bai, fastq.star_bam, "samtools index " + fastq.star_bam)

    # Let's run some rseqc
    ref_seq_bed = join(genome_annotation_dir, "hg19_RefSeq.bed")
    fastq.read_dist_report = join(rseqc_dir, fastq.fastq_sample, "read_distribution.txt")
    dm.add(
        fastq.read_dist_report,
        [fastq.star_bam, ref_seq_bed],
        "{} && read_distribution.py -i {} -r {} > {}".format(activate_str, fastq.star_bam,
                                                             ref_seq_bed, fastq.read_dist_report)
    )

# QC summary
multiqc_report = join(results_dir, "multiqc", "multiqc_report.html")
dm.add(multiqc_report, chain(*[(f.fastqc_report, f.star_bam, f.read_dist_report) for f in fastqs]),
       "multiqc -d -dd 1 -f --outdir {} {} {} {}".format(dirname(multiqc_report), fastqc_dir,
                                                         alignment_dir,
                                                         rseqc_dir))

# combine experimental samples
exp_samples = collate_fastqs_to_experimental_samples(fastqs)
for sample in exp_samples:
    star_bams = [f.star_bam for f in sample.fastqs]
    sample.star_bam_merged = join(alignment_dir, "star_by_experimental_sample", sample.name,
                                  sample.name + ".bam")
    merge_dir = dirname(sample.star_bam_merged)
    dm.add(
        sample.star_bam_merged, star_bams, ''.join([
            "cd {}".format(merge_dir),
            " && rm -rf {} && ".format(sample.star_bam_merged),
            " && ".join(["ln -f -s {} {}".format(f.star_bam, f.fastq_sample + ".bam") for f in sample.fastqs]),
            " && ",
            "samtools merge -r {} *.bam".format(sample.star_bam_merged),
            " && samtools index " + sample.star_bam_merged
        ])
    )

    # HT-seq - count reads on exons
    # sample.htseq_counts = join(htseq_dir, sample.name, sample.name + ".counts")
    # dm.add(
    #     sample.htseq_counts,
    #     [sample.star_bam_merged, ensembl_annotations_with_chr_no_zip],
    #     activate_str +
    #     " && htseq-count --format bam --stranded no --quiet"
    #     " {} {} > {}".format(sample.star_bam_merged,
    #                          ensembl_annotations_with_chr_no_zip,
    #                          sample.htseq_counts)
    # )

    # try featureCounts - count reads on exons
    sample.featureCounts_counts = join(featureCounts_dir, sample.name, "exons", sample.name + ".exon.counts")
    dm.add(
        sample.featureCounts_counts,
        [sample.star_bam_merged, ensembl_annotations_with_chr_no_zip],
        " featureCounts"
        " -t exon -g gene_id -a {} -o {} {}".format(ensembl_annotations_with_chr_no_zip,
                                                    sample.featureCounts_counts,
                                                    sample.star_bam_merged, )
    )

    # reformat counts to simpler format
    sample.featureCounts_counts_simple = sample.featureCounts_counts + ".simple_counts"
    dm.add(
        sample.featureCounts_counts_simple,
        sample.featureCounts_counts,
        "grep -v '^#' {} | cut -f 1,7 > {}".format(sample.featureCounts_counts, sample.featureCounts_counts_simple)
    )

sample_multiqc_report = join(results_dir, "multiqc_per_experimental_sample", "multiqc_report.html")
dm.add(sample_multiqc_report,
       chain(*[(s.featureCounts_counts,) for s in exp_samples]),
       "multiqc -f --outdir {} {} {}".format(dirname(sample_multiqc_report), htseq_dir,
                                             featureCounts_dir))

# collate counts into one csv
exp_samples = sorted(exp_samples, key=lambda x: x.name)
fc_counts = [s.featureCounts_counts_simple for s in exp_samples]
fc_merged_counts = join(featureCounts_dir, "all_sample_counts.csv")
assert len(fc_counts) >= 2
dm.add(
    fc_merged_counts, fc_counts,
    "rm -f {}".format(fc_merged_counts) +
    " && join {} {}".format(*fc_counts[0:2]) + "".join(
        [" | join - " + f for f in fc_counts[2:]]) + " >> " + fc_merged_counts
)

# run R analysis
flag = join(htseq_dir, "analysis.R.flag")
# dm.add(flag, htseq_merged_counts, "Rscript analysis.R")

dm.execute()

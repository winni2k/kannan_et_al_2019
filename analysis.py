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


class DMWorker(object):
    def __init__(self, dm):
        self.dm = dm

    def unzip_to(self, zipped, unzipped):
        self.dm.add(unzipped, zipped,
                    "gzip -dc {} > {}".format(zipped, unzipped))


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


MAX_READ_LENGTH = 75
MAPPING_THREADS = 4

args = get_dm_arg_parser('RNA-seq analysis').parse_args()
dm = DistributedMake(args_object=args, keep_going=True)
dmw = DMWorker(dm)

activate_str = ". activate rna_seq_py27"

# This is the location where Kiran downloaded the FASTQ files to
rna_seq_dir = abspath("data/AdHa160517-40613574")
genome_annotation_dir = abspath("data/QC/annotation")
ensembl_annotations = abspath("data/annotations/ensembl/Homo_sapiens.GRCh38.89.gtf.gz")
ensembl_reference = abspath("data/reference/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
refseq_reference = abspath("data/reference/refseq/GCF_000001405.36_GRCh38.p10_genomic.fna.gz")
refseq_annotations = abspath("data/annotations/refseq/GCF_000001405.36_GRCh38.p10_genomic.gff.gz")
go_terms = abspath("data/gene_sets/Human_GOALL_with_GO_iea_June_01_2017_symbol.gmt")

results_dir = abspath("results")
fastqc_dir = join(results_dir, "fastqc")
rseqc_dir = join(results_dir, "rseqc")
alignment_dir = join(results_dir, "alignment")
htseq_dir = join(results_dir, "htseq")
featureCounts_dir = join(results_dir, "featureCounts")
ensembl_annot_dir = join(results_dir, "annotations", "ensembl")
ensembl_ref_dir = join(results_dir, "reference", "ensembl")
gene_set_dir = join(results_dir, "gene_sets")

ensembl_reference_unzip = join(ensembl_ref_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
dmw.unzip_to(ensembl_reference, ensembl_reference_unzip)

ensembl_annotations_unzip = join(ensembl_annot_dir, "Homo_sapiens.GRCh38.89.gtf")
dm.add(ensembl_annotations_unzip, ensembl_annotations,
       "gzip -dc {} > {}".format(ensembl_annotations, ensembl_annotations_unzip))

ensembl_annotations_real_genes = join(ensembl_annot_dir,
                                      "Homo_sapiens.GRCh38.89.protein_coding.gtf")
dm.add(ensembl_annotations_real_genes, ensembl_annotations_unzip,
       (' grep -e \'^#\' -e \'gene_biotype "protein_coding"\' {}'
        ' > {}').format(ensembl_annotations_unzip, ensembl_annotations_real_genes)
       )

# create a new STAR index from GRCh37
star_ref_flag = join(results_dir, "reference/GRCh38.89_star/flag")
star_ref_dir = dirname(star_ref_flag)
dm.add(star_ref_flag, [ensembl_reference_unzip, ensembl_annotations],
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
fastqs = gather_fastq_info(rna_seq_dir)
for fastq in fastqs:
    fastq.fastqc_report = join(fastqc_dir,
                               fastq.experimental_sample,
                               basename(fastq.fastq).rstrip(".fastq.gz") + "_fastqc.html")
    out_dir = dirname(fastq.fastqc_report)
    dm.add(fastq.fastqc_report, fastq.fastq,
           "fastqc --outdir {} {}".format(out_dir, fastq.fastq))

# Let's map the reads using hisat2 - oh wait hisat2 quality sucks
# STAR it is.
for fastq in fastqs:
    fastq.star_bam = join(alignment_dir, "star", fastq.fastq_sample,
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

# QC summary
multiqc_report = join(results_dir, "multiqc", "multiqc_report.html")
dm.add(multiqc_report, chain(*[(f.fastqc_report, f.star_bam) for f in fastqs]),
       "multiqc -d -dd 1 -f --outdir {} {} {}".format(dirname(multiqc_report), fastqc_dir,
                                                      alignment_dir))

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
            " && ".join(["ln -f -s {} {}".format(f.star_bam, f.fastq_sample + ".bam") for f in
                         sample.fastqs]),
            " && ",
            "samtools merge -r {} *.bam".format(sample.star_bam_merged),
            " && samtools index " + sample.star_bam_merged
        ])
    )

    # try featureCounts - count reads on exons
    sample.featureCounts_counts = join(featureCounts_dir, sample.name, "exons",
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

sample_multiqc_report = join(results_dir, "multiqc_per_experimental_sample", "multiqc_report.html")
dm.add(sample_multiqc_report,
       chain(*[(s.featureCounts_counts,) for s in exp_samples]),
       "multiqc -f --outdir {} {} {}".format(dirname(sample_multiqc_report), htseq_dir,
                                             featureCounts_dir))

# calculate gene lengths
fc_gene_lengths = join(featureCounts_dir, "gene_lengths.tsv")
dm.add(fc_gene_lengths, exp_samples[0].featureCounts_counts,
       ("grep '^ENSG' {}"
        " | cut -f1,3,4"
        " | perl -ane 'print $$F[0] . \"\\t\" . ($$F[2] - $$F[1]) . \"\\n\"' > {} ").format(
           exp_samples[0].featureCounts_counts, fc_gene_lengths))

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

# convert go terms to go id only
go_terms_only_symbols = join(gene_set_dir,
                             "Human_GOALL_with_GO_iea_June_01_2017_symbol.gmt")
dm.add(go_terms_only_symbols, go_terms,
       "perl -pe 's/.*\%//' < {} > {}".format(go_terms, go_terms_only_symbols))

# run R analysis
flag = join(htseq_dir, "analysis.R.flag")
# dm.add(flag, htseq_merged_counts, "Rscript analysis.R")

dm.execute()

from os.path import join, basename
from shlex import quote
from urllib.parse import urlparse

from dmpy.distributedmake import DistributedMake, get_dm_arg_parser


class Downloader(object):
    def __init__(self, dm):
        self.dm = dm

    def get_dir(self, url, out_dir):
        self.dm.add(
            join(out_dir, "flag"),
            None,
            ('cd {} '
             '&& wget -r -nH -nd -np -R index.html* {} '
             '&& touch flag').format(quote(out_dir), url)
        )

    def get_file(self, url, output):
        command = "curl {} -o {}".format(quote(url), quote(output))
        self.dm.add(output, None, command)

    def get_file_to_dir(self, url, output_dir):
        output = join(output_dir, basename(urlparse(url).path))
        self.get_file(url, output)


args = get_dm_arg_parser('RNA-seq analysis').parse_args()
dm = DistributedMake(args_object=args)
down = Downloader(dm)

# down.get_dir(
#     "https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/reference/hg19_Gencode14.overhang75/",
#     "data/reference/hg19_Gencode14.overhangt75_star"
# )
# down.get_dir(
#     "https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/QC/annotation/",
#     "data/QC/annotation")

down.get_file_to_dir(
    'ftp://ftp.ensembl.org/pub/release-89/gtf/homo_sapiens/Homo_sapiens.GRCh38.89.gtf.gz',
    "data/annotations/ensembl"
)

down.get_file_to_dir(
    'ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz',
    'data/reference/ensembl'
)

down.get_file_to_dir(
    'http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/GO/Human_GOALL_with_GO_iea_June_01_2017_symbol.gmt',
    'data/gene_sets'
)

broad_bundle_urls = [
    'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz',
    'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi',
    'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
    'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi',
    'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz',
    'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi',
]
for url in broad_bundle_urls:
    down.get_file_to_dir(url, 'data/broad/bundle/hg38')

dm.execute()

# 2017-07-01 downloaded from email from pavitra: data/sample_meta_data.tsv

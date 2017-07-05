from os.path import join
from shlex import quote

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


args = get_dm_arg_parser('RNA-seq analysis').parse_args()
dm = DistributedMake(args_object=args)
down = Downloader(dm)

down.get_dir(
    "https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/reference/hg19_Gencode14.overhang75/",
    "data/reference/hg19_Gencode14.overhangt75_star"
)
down.get_dir(
    "https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/QC/annotation/",
    "data/QC/annotation")

down.get_file(
    'ftp://ftp.ensembl.org/pub/release-89/gtf/homo_sapiens/Homo_sapiens.GRCh38.89.gtf.gz',
    "data/annotations/ensembl/Homo_sapiens.GRCh37.89.gtf.gz"
)

down.get_file(
    'ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz',
    'data/reference/ensembl/Homo_sapiens.GRCh38.dna.toplevel.fa.gz'
)
# 2017-07-01 downloaded from email from pavitra: data/sample_meta_data.tsv


dm.execute()

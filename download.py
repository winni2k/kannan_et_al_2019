from os.path import join
from shlex import quote

from dmpy.distributedmake import DistributedMake, get_dm_arg_parser


def get_dir(dm, url, out_dir):
    dm.add(
        join(out_dir, "flag"),
        None,
        ('cd {} '
         '&& wget -r -nH -nd -np -R index.html* {} '
         '&& touch flag').format(quote(out_dir), url)
    )


def get_file(dm, url, output):
    command = "curl {} -o {}".format(quote(url), quote(output))
    dm.add(output, None, command)


args = get_dm_arg_parser('RNA-seq analysis').parse_args()
dm = DistributedMake(args_object=args)

get_dir(dm,
        "https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/reference/hg19_Gencode14.overhang75/",
        "data/reference/hg19_Gencode14.overhangt75_star")
get_dir(dm,
        "https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/QC/annotation/",
        "data/QC/annotation")

get_file(dm,
         'ftp://ftp.ensembl.org/pub/grch37/release-89/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz',
         "data/annotations/ensembl/Homo_sapiens.GRCh37.87.gtf.gz")

dm.execute()

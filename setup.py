from setuptools import setup

setup(
    name='genpipes',
    version='4',
    packages=['bfx', 'core', 'utils', 'pipelines', 'pipelines.epiqc', 'pipelines.covseq', 'pipelines.dnaseq', 'pipelines.hicseq', 'pipelines.rnaseq', 'pipelines.chipseq', 'pipelines.nanopore', 'pipelines.methylseq', 'pipelines.tumor_pair', 'pipelines.ampliconseq', 'pipelines.rnaseq_light', 'pipelines.nanopore_covseq', 'pipelines.dnaseq_high_coverage', 'pipelines.rnaseq_denovo_assembly'],
    package_dir={'': 'genpipes'},
    url='https://bitbucket.org/mugqic/genpipes/',
    license='GNU Lesser General Public License',
    author='Canadian Center for Comoutational Genomics',
    author_email='pipelines@computationalgenomics.ca.',
    description='Several bioinformatics pipelines developed at McGill University Genome Centre'
)

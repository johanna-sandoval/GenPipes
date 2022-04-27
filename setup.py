import setuptools

setuptools.setup(
    name='genpipes',
    version='4',
    packages=setuptools.find_packages(),
    url='https://bitbucket.org/mugqic/genpipes/',
    license='GNU Lesser General Public License',
    author='Canadian Center for Comoutational Genomics',
    author_email='pipelines@computationalgenomics.ca.',
    description='Several bioinformatics pipelines developed at McGill University Genome Centre',
    entry_points = {
    "console_scripts": [
        "genpipes = genpipes.__main__:main"] },
    install_requires = ["packaging >=   20.9"]
)

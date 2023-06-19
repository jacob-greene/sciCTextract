from pathlib import Path

from setuptools import setup, find_packages, find_namespace_packages

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

entry_points = {
    'console_scripts': [
        'sciCTextract = sciCTextract.extract:main'
    ]
}

install_requires = [
    "biopython>=1.79",
]

setup(
    name='sciCTextract',
    version='0.2',
    author='Matt Fitzggibon',
    license='MIT',
    author_email='mfitzgib@fredhutch.org',
    description='Simple demultiplexer for sciCUT&Tag data',
    include_package_data=True,
    packages=find_namespace_packages(),
    install_requires=install_requires,
    url='https://github.com/mfitzgib/sciCTextract',
    entry_points=entry_points,
    long_description=long_description,
    long_description_content_type="text/markdown"
)

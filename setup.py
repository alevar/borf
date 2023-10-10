from setuptools import setup, Extension

extensions = [
    Extension("borf.classes.treader", ["borf/classes/treader.pyx"]),
    Extension("borf.classes.transcript", ["borf/classes/transcript.pyx"]),
    Extension("borf.classes.txgroup", ["borf/classes/txgroup.pyx"]),
]

setup(
    name="borf",
    version="0.0.1",
    author="Ales Varabyou, Zayn Zaidi",
    author_email="ales.varabyou@jhu.edu",
    description="Selecting most representative ORF by clustering via ILPI",
    url="https://github.com/alevar/borf",
    install_requires=["argparse","unittest"],
    python_requires='>=3.6',
    packages=['borf'],
    entry_points={'console_scripts': ['borf = borf.main:main'], },
    package_dir={"borf":"borf"},
    ext_modules=cythonize(extensions),
)

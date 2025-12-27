"""Setup script for Lipid-Protein Interaction Analyzer"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="lipid-protein-analyzer",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A comprehensive tool for analyzing lipid-protein interactions from MD simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/lipid-protein-analyzer",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19.0",
        "pandas>=1.1.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "MDAnalysis>=2.0.0",
        "tqdm>=4.50.0",
        "scipy>=1.5.0",
    ],
    entry_points={
        "console_scripts": [
            "lipid-analyzer-stage1=stage1_contact_analysis.main:main",
        ],
    },
)
#!/usr/bin/env python
"""Setup script for LIPAC package"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="lipac",
    version="1.0.0",
    author="Takeshi Sato",
    author_email="takeshi.sato@mb.kyoto-phu.ac.jp",
    description="Lipid-Protein Analysis with Causal Inference for MD simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/takeshi-sato-dev/lipac",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "MDAnalysis>=2.0.0",
        "pymc>=5.0.0",
        "arviz>=0.14.0",
        "tqdm>=4.60.0",
        "scikit-learn>=1.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
            "sphinx",
            "sphinx-rtd-theme",
        ],
    },
    entry_points={
        "console_scripts": [
            "lipac-stage1=stage1_contact_analysis.main:main",
            "lipac-stage2=stage2_contact_analysis.main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.md", "*.txt", "*.yml", "*.yaml"],
    },
    zip_safe=False,
)

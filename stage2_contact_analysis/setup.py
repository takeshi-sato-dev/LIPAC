"""
Setup script for BayesianLipidAnalysis package
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="bayesian_lipid_analysis",
    version="1.0.0",
    author="Takeshi Sato",
    author_email="your.email@example.com",  # Update with actual email
    description="Bayesian analysis of lipid-protein interactions from molecular dynamics simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/bayesian_lipid_analysis",  # Update with actual URL
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "bayesian_lipid_analysis=bayesian_lipid_analysis_main:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    keywords="molecular-dynamics bayesian-statistics lipid-protein-interactions membrane-proteins",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/bayesian_lipid_analysis/issues",
        "Source": "https://github.com/yourusername/bayesian_lipid_analysis",
        "Documentation": "https://bayesian_lipid_analysis.readthedocs.io",
    },
)
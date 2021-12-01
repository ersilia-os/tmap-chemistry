from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf8") as fh:
    long_description = fh.read()

install_requires = [
    "click"
]


setup(
    name="tmapchem",
    version="0.0.1",
    author="Miquel Duran-Frigola",
    author_email="miquel@ersilia.io",
    url="https://github.com/ersilia-os/tmap-chemistry",
    description="Fast visualization of the chemical space",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.7",
    install_requires=install_requires,
    entry_points={"console_scripts": ["tmapchem=tmapchem.cli:cli"]},
    packages=find_packages(exclude=("utilities")),
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Data Visualization",
    ],
    keywords="chemical-space visualisation",
    project_urls={"Source Code": "https://github.com/ersilia-os/tmap-chem",},
    include_package_data=True,
)
# Visualize the chemical space with TMAP

This Python library is a simple wrapper for the fantastic [TMAP library](https://github.com/reymond-group/tmap) developed by the [Reymond Group](https://www.gdb.unibe.ch/). All credit goes to the original authors: here we just provide a quick interface with limited functionality. If you want to become a frequent user of TMAP, please refer the [documenation of the TMAP library](https://tmap.gdb.tools/), which includes many and varied examples. They even provide a [web server](https://try-tmap.gdb.tools/) to try the tool for up to 15k molecules.

With the current library, we want to accomplish two goals:
1. Connect TMAP to the [Ersilia Model Hub](https://ersilia.io). This means annotating chemical libraries with any of the properties available from our hub.
2. Merge small libraries (a few thousand compounds) with larger libraries like ChEMBL or Natural Products (e.g. Coconut), so that a unified view is accomplished.

## Installation

We recommend that you create a conda environment:

```bash
# create conda environment
conda create -n tmapchem python=3.7

# activate
conda activate tmapchem
```

Then, simply clone this repository and install it using pip:

```bash
# clone repository
git clone git@github.com:ersilia-os/tmap-chemistry.git

# install with pip
cd tmap-chemistry
python -m pip install -e .
```

The main dependencies are, obviously, [tmap](https://github.com/reymond-group/tmap) and [ersilia](https://github.com/ersilia-os/ersilia).

```bash
# install tmap and faerun
conda install -c tmap tmap
python -m pip install faerun

# install ersilia
python -m pip install git+ssh://git@github.com/ersilia-os/ersilia.git
```

## Usage

You can run an example as follows:

```bash
tmapchem -n my_project -i input_file.csv -p params.yaml -o output_folder
```

This command will produce an output folder containing an HTML file. Simply double click and the chemical space will be displayed in your folder!

For more help, type:

```bash
tmapchem --help
```

### The parameters file

You can specify a parameters file (named `params.yaml`in the example below) containing model identifiers or slugs from the Ersilia Model Hub. The file format should be like this:
```yaml
ersilia-hub: 
 - chemistry-default-descriptors
```

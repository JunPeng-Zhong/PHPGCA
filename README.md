# Overview

## PHPGCA

PHPGCA is a phage host prediction tool that employs graph convolution and contrastive learning.

## Required Dependencies

### Clone the Repository

To get started, clone the repository using the following command:

```bash
git clone https://github.com/JunPeng-Zhong/PHPGCA.git
```

### Setting Up the Environment

You'll need to set up the environment using either conda or miniconda. Navigate to the project directory and create the environment with the required dependencies by running:

```bash
cd PHPGCA
conda env create -f requirements_conda.yaml -n PHPGCA
```

Once the environment is created, activate it using:

```bash
conda activate PHPGCA
```

*Note:* The tool is based on `pyg`, if the cuda version is uncompatity, please follow the installation step in the [Installation — pytorch_geometric documentation](https://pytorch-geometric.readthedocs.io/en/latest/install/installation.html)

Due to GitHub limitations, the datasets are hosted on [Google Drive link](https://drive.google.com/drive/folders/1rMcOsaZPozy3f0USxqkH5Orhn0KfFjcX?usp=sharing) . Firstly, download the `blast_db1/` and move it to `PHPGCA/` 

You should download the datasets and then move them to the `datas/` directory.

As an example, we'll consider the CHERRY dataset:

```bash
python switch_dataset.py --dataset CHERRY
```

## Usage

### Switching Datasets

You can easily switch between datasets by using the `switch_dataset.py` script. Download the desired dataset, move it to the `datas/` directory, and then execute:

```bash
python switch_dataset.py --dataset [dataset_name]
```

### Predicting Hosts for Viruses

To reproduce the results, use the provided scripts for top 1 and top 5 predictions:

For top 1 prediction:

```bash
./pretrain.sh
```

For top 5 predictions:

```bash
./pretrain_5.sh
```

**Example**

```bash
python switch_dataset.py --dataset CHERRY
./pretrain.sh
```

If you want to test the host of a custom virus, ensure you have a fasta file named `test_contigs.fa` containing viral sequences. Then, use the following commands:

```bash
./construct_graph.sh
./retrain.sh
```

**Output**

The predictions for test viruses will be saved in the `final_prediction.csv` file, with the `contig` column representing the input accession.

### Adding Additional Training Data

To include your own data for **training**, follow these steps:

- Append the virus sequences to `dataset/nucl.fasta`.
- Add the corresponding host taxonomy to `dataset/virus.csv`, including `Accession` and `Species`.
- Place prokaryotic genomes in the `prokaryote/` directory and update `dataset/prokaryote.csv` with their taxonomy details (`Accession` and `Species`).

This will allow you to extend the training data with your custom information.

## Contact

Feel free to ask if you have any further questions: g597234159@gmail.com
Welcome to the ligand binding site predictor.

Before running the predictor algorithm, install all the necessary requirements by typing:

``` pip install -r requirements.txt```

To run the predictor simply type the following command in the command line:

```python3 main.py path/to/pdb/file```

The machine learning predictor based on a random forest algorithm will then train itself on its data and output names of residues that are possibly ligand binding for the input file provided. The **chimera_cmds.cmd** file is created as well, which we will use in the next step.

Next, open the input file in Chimera and do the following steps:

- open command line interface in chimera

- type:

``` open /path/to/the/created/chimera_cmds.cmd```

It will select the residues that are predicted to be ligand binding.



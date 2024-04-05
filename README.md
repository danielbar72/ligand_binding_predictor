Welcome to the ligand binding site predictor.

Before running the predictor algorithm, install all the necessary requirements by typing:

``` pip install -r requirements.txt```

To run the predictor simply type the following command in the command line:

```python3 main.py path/to/pdb/file```

The machine learning predictor based on a random forest algorithm will then train itself on its data and output names of residues that are possibly ligand binding for the input file provided. 



For testing purposes, run the following command:

```python3 main.py test/3rt6.pdb```



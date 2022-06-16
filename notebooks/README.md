# Tutorial on GCP & Python bits
This requires some libraries to be installed. These are at the root level of the repository in a requirements.txt file

Assuming you have conda installed, good practice would be to create a virtual environment per repository

```bash
$ conda create -n rd_venv python=3.10
$ conda activate rd_venv
# this should show the active environment on the prompt
$ (rd_venv) conda install --file requirements.txt
```

After the package installation completes, you should be able to run the accompanying notebook with

```bash
jupyter notebook notebooks/get_cram_paths.ipynb
```
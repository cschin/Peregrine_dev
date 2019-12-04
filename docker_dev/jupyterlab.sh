#!/bin/bash
. /opt/conda/etc/profile.d/conda.sh
conda activate peregrine
jupyter-lab --ip="*" --allow-root --no-browser --port 8888 --NotebookApp.token='' --NotebookApp.disable_check_xsrf=True /wd

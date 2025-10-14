set -e
jupyter nbconvert --to notebook --inplace --execute block_aligner_bench_vis.ipynb --allow-errors
jupyter trust block_aligner_bench_vis.ipynb
jupyter nbconvert --to notebook --inplace --execute block_aligner_accuracy_vis.ipynb --allow-errors
jupyter trust block_aligner_accuracy_vis.ipynb

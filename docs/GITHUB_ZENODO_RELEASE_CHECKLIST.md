# GitHub And Zenodo Release Checklist

## Before Pushing To GitHub

- Confirm the repository contains code and documentation, not large downloaded data.
- Run `python scripts/build_code_release_manifest.py`.
- Check that `docs/code_file_manifest.tsv` and `docs/script_manifest.tsv` were updated.
- Run `python -m compileall scripts` if the local environment has the required Python version.
- Confirm `.gitignore` excludes `data/`, `results/`, `figures/`, `manuscript/Final figures/`, manuscript drafts, and caches.
- Confirm the license choice. The current prepared code license is MIT.
- Update `CITATION.cff` and `.zenodo.json` with the final GitHub repository URL after the repository exists.

## Suggested GitHub Files

- `README.md`
- `LICENSE`
- `CITATION.cff`
- `.zenodo.json`
- `CODE_AVAILABILITY.md`
- `requirements.txt`
- `environment.yml`
- `config/`
- `scripts/`
- `docs/`

## Suggested Zenodo Workflow

1. Push the clean code repository to GitHub.
2. Connect the GitHub repository to Zenodo.
3. Create a GitHub release, for example `v0.1.0`.
4. Let Zenodo archive the GitHub release automatically.
5. Add the Zenodo DOI to `CITATION.cff`, `README.md`, and the manuscript code-availability section after Zenodo assigns it.

## Do Not Deposit In The Code Archive

- Raw or processed controlled-access dbGaP matrices.
- Full local `data/` cache.
- Full local `results/` cache unless a specific journal requires it separately.
- Manuscript `.docx` drafts.
- Local browser/cache folders.

## Optional Separate Deposit

If desired, final figure JPEGs and supplementary tables can be deposited as a separate Zenodo "dataset" record, distinct from the code/software record. That record should use data/table licensing rather than the MIT code license.

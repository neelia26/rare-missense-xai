# Minimal Git Workflow

The repository should preserve code and provenance without publishing the large source databases.

```bash
git init
git add README.md SOURCE_FILES_NEEDED.md DATA_DICTIONARY.md GIT_WORKFLOW.md RUN_LOG_TEMPLATE.md
git add VERSION requirements.txt run_pipeline.py config scripts src tests .gitignore
git commit -m "Add governed rare missense dataset pipeline"
```

After a successful stage, review its audit before committing code or small audit records:

```bash
python run_pipeline.py status
git status
```

Do not use `git add .` until you understand every untracked file. Raw databases, API caches and generated patient/variant-scale tables should remain outside the public repository. For each external source, preserve its release, download URL, retrieval date and SHA-256 in a source manifest or run log.

Recommended branches:

- `main`: reviewed scripts and documentation;
- `pilot`: small-cohort testing changes; and
- short feature branches for individual enrichment or audit stages.

Tag the exact code used for the final reconstruction:

```bash
git tag -a dataset-pipeline-v1.0.0 -m "Pipeline used for governed dataset rebuild"
```


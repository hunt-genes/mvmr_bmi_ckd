# Uploading this to GitHub

The repository `hunt-genes/mvmr_bmi_ckd` already exists (currently empty).
From your local computer, after downloading and unzipping this folder:

```bash
cd mvmr_bmi_ckd
git init
git add .
git commit -m "Analysis code for tissue-partitioned MR of BMI on kidney function"
git branch -M main
git remote add origin https://github.com/hunt-genes/mvmr_bmi_ckd.git
git push -u origin main
```

`data/raw/` is excluded by `.gitignore`, so your input GWAS files stay off
GitHub. Only the code, README, and empty folder placeholders are pushed.

If the remote already has a commit (e.g. an auto-created README), run
`git pull --rebase origin main` before `git push`.

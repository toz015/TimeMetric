# History Purge Plan — borrowed implementation blobs

**Date prepared:** 2026-09-24
**Status:** **PREPARED ONLY. NOT EXECUTED.** No rewrite, push, force-push, tag,
release, or CRAN submission has been performed.
**Requires:** explicit maintainer approval before any command below is run.
**Verified with:** `git filter-repo --dry-run` in a throwaway mirror clone, which
writes export streams and rewrites nothing. `git-filter-repo` version
`a40bce548d2c`, at `/opt/homebrew/bin/git-filter-repo`.

---

## 1. Purpose and scope

Deleting `R/pam.rsph.R` and `R/pam.schemper.R` in commit `644da37` stops the
package *distributing* borrowed code, which is what CRAN and JOSS assess. It does
not remove that code from git history, where it remains reachable.

This purge removes **borrowed implementation code and borrowed generated
documentation** from all reachable history. It does **not** remove discussion:
the strings `r_e`, `r_sh`, `rsph` and `schemper` in commit messages, NEWS, the
findings log, the audit and the design documents are legitimate development
history and are preserved deliberately.

## 2. Regenerated purge-path list

Derived by scanning the **content** of every blob in every commit for
implementation-specific fragments, restricted to `R/` and `man/`. It was not
copied from the design spec — and regenerating it mattered: **the spec listed 9
paths and the correct number is 10.** `man/pam.rsph_metric.Rd` was missing.

Scan command:

```bash
git rev-list --all | while read c; do
  git grep -l -E 'tempi\.eventi|num\.sogg|ind\.censura|stima\.surv|f\.assegna\.surv|newvare|incom\.sum|n\.incom|Re\.imp|Re\.fix|r2nw|dRti|sen0' "$c" -- R man 2>/dev/null
done | awk -F: '{print $2}' | sort | uniq -c | sort -rn
```

| # | Path | Blobs | Why |
|---|---|---|---|
| 1 | `R/pam.rsph.R` | 11 | 92% verbatim from the unlicensed `Re.r` |
| 2 | `R/pam.re.R` | 4 | Same code under its original name; `pam.rsph.R` was renamed from it at 97% similarity |
| 3 | `R/pam.rsph.metric.R` | 1 | Duplicate `R_E` implementation |
| 4 | `R/pam.rsph_metric.R` | 4 | **Distinct path** from #3 — underscore, not dot |
| 5 | `R/pam.schemper.R` | 6 | 91% verbatim from GPL-2 `survAUC::schemper()` |
| 6 | `man/pam.rsph.Rd` | 1 | Roxygen output reproducing the borrowed signature |
| 7 | `man/pam.re.Rd` | 3 | As above |
| 8 | `man/pam.rsph.metric.Rd` | 1 | As above |
| 9 | `man/pam.rsph_metric.Rd` | 1 | **Missing from the spec.** Documents `pam.rsph_metric(predicted_data, survival_time, status, start_time)` |
| 10 | `man/pam.schemper.Rd` | 2 | Documents `pam.schemper(train.fit, traindata, newdata)` |

**34 distinct blobs** across the 10 paths.

The `man/*.Rd` pages are included because they are roxygen-generated *from* the
borrowed sources and reproduce their signatures and argument documentation,
which is copied documentation.

### Deliberately NOT purged

| Path | Reason |
|---|---|
| `docs/superpowers/r-e-implementation-audit.md` | Quotes the borrowed identifiers as diagnostic evidence |
| `docs/superpowers/findings.md` | Records how the licensing problem was found |
| `docs/superpowers/specs/2026-09-24-…-design.md`, `plans/2026-09-24-…md` | Quote fragments to define the verification |
| `tests/testthat/test-rsph-dispatch.R`, `test-r-e-reference.R`, `test-r-sh-definition.R`, `_snaps/rsph-dispatch.md` | Original work. Zero matches for any `survAUC`-port identifier; they name `dRti` and `my.survfit` only as assertion targets, e.g. `expect_identical(names(res), c("times","Rti","dRti"))`. The snapshot holds base64-serialised *output*, not source |
| `R/tm_survival_eval.R`, `R/tm_fit_and_eval.R`, `R/pam.predicted_survial_eval.R` and similar | These merely **called** the borrowed functions. A broad `rsph|schemper` substring scan matches them; purging them would destroy legitimate package history. This is why the narrow fragment set is the correct discriminator |

## 3. Refs and blast radius

```
refs/heads/joss-revision    1f01dd2   local only, never pushed
refs/heads/main             1f54e6e   PUBLISHED - identical to origin/main
refs/remotes/origin/main    1f54e6e
refs/remotes/origin/HEAD    1f54e6e
tags                        (none)
```

* **106 commits** on all refs; **19** touch the purge paths. Those 19 and every
  descendant get new hashes, so in practice **all 106 commits are rewritten**.
* `git filter-repo` rewrites **all refs**, not a selected branch.

### ⚠️ The consequential risk

**`R/pam.rsph.R` and `R/pam.schemper.R` are present on `origin/main`**, and local
`main` is identical to `origin/main`. So the purge cannot be confined to the
unpushed branch:

* `main` must be rewritten and **force-pushed**, rewriting published history on
  `github.com/toz015/TimeMetric`.
* Every existing clone, fork and open PR becomes incompatible. Collaborators
  must re-clone; they cannot merge or rebase across the rewrite.
* Anyone who already cloned retains the borrowed code locally. The purge limits
  future distribution; it cannot retract what has already been fetched.
* GitHub may retain unreferenced objects until its own GC runs. Blobs can stay
  reachable via the API by SHA for a period. **GitHub Support must be asked to
  purge cached views**, and any fork must be deleted, or the blobs survive there.

This is a decision about a published repository, not a local cleanup, and is why
nothing here runs without explicit approval.

## 4. Backup and rollback

Take **both** backups. They are cheap and the operation is otherwise
irreversible.

```bash
# 4.1 Full mirror backup, outside the working repo
git clone --mirror https://github.com/toz015/TimeMetric.git \
  ~/timemetric-backup-$(date +%Y%m%d)-remote.git

cd /Users/wanghd/Documents/workplace/TimeMetrics/TimeMetric
git bundle create ~/timemetric-backup-$(date +%Y%m%d)-local.bundle --all

# 4.2 Record every ref so any branch can be restored by hash
git for-each-ref --format='%(objectname) %(refname)' \
  > ~/timemetric-refs-$(date +%Y%m%d).txt

# 4.3 Record the blobs that must become unreachable
for p in R/pam.rsph.R R/pam.re.R R/pam.rsph.metric.R R/pam.rsph_metric.R \
         R/pam.schemper.R man/pam.rsph.Rd man/pam.re.Rd \
         man/pam.rsph.metric.Rd man/pam.rsph_metric.Rd man/pam.schemper.Rd; do
  for c in $(git log --all --pretty='%H' -- "$p"); do
    b=$(git rev-parse --quiet --verify "$c:$p" 2>/dev/null) && echo "$p $b"
  done
done | sort -u > ~/timemetric-purge-blobs-$(date +%Y%m%d).txt
wc -l ~/timemetric-purge-blobs-*.txt    # expect 34
```

**Rollback before any push** — local only, fully recoverable:

```bash
# filter-repo leaves the pre-rewrite refs under refs/original/
git for-each-ref --format='%(refname)' refs/original/
git update-ref refs/heads/main $(git rev-parse refs/original/refs/heads/main)
git update-ref refs/heads/joss-revision $(git rev-parse refs/original/refs/heads/joss-revision)
# or discard the rewritten repo entirely and restore:
git clone ~/timemetric-backup-YYYYMMDD-local.bundle TimeMetric-restored
```

**Rollback after a force-push** — only from the mirror backup:

```bash
cd ~/timemetric-backup-YYYYMMDD-remote.git
git push --mirror https://github.com/toz015/TimeMetric.git
```

Collaborators must re-clone again. Treat post-push rollback as damage control,
not a routine undo.

## 5. Pre-rewrite verification

```bash
# clean tree, known commit
git status --short                      # must be empty
git rev-parse HEAD                      # record it

# full suite green before touching history
Rscript -e 'testthat::test_local()'     # expect 498 pass / 0 fail / 0 skip

# the 34 blobs are currently reachable (sanity: the purge has something to do)
while read p b; do
  git cat-file -e "$b" 2>/dev/null || echo "ALREADY UNREACHABLE: $p $b"
done < ~/timemetric-purge-blobs-YYYYMMDD.txt
```

## 6. The rewrite — NOT EXECUTED

Run in a **fresh mirror clone**, which is what `filter-repo` expects; do not run
it in the working repository.

```bash
git clone --mirror https://github.com/toz015/TimeMetric.git TimeMetric-rewrite.git
cd TimeMetric-rewrite.git

git filter-repo \
  --invert-paths \
  --path R/pam.rsph.R \
  --path R/pam.re.R \
  --path R/pam.rsph.metric.R \
  --path R/pam.rsph_metric.R \
  --path R/pam.schemper.R \
  --path man/pam.rsph.Rd \
  --path man/pam.re.Rd \
  --path man/pam.rsph.metric.Rd \
  --path man/pam.rsph_metric.Rd \
  --path man/pam.schemper.Rd
```

Add `--force` only if filter-repo refuses because the clone is not considered
fresh. Note that `filter-repo` **removes the `origin` remote** after rewriting,
by design, so the push in §8 re-adds it explicitly.

**This should be merged with the separately planned key-material purge** (the
committed `q` / `q.pub` ed25519 pair, spec step 8) so history is rewritten
**once**, not twice. Add their paths to the same invocation.

### Dry-run result (this was run; it rewrites nothing)

```
Parsed 106 commits
New history written in 0.04 seconds; now repacking/cleaning...
NOTE: Not running fast-import or cleaning up; --dry-run passed.
```

Comparing the two streams, **all 10 paths disappear as file operations**:

| Path | file ops before | after |
|---|---|---|
| `pam.rsph.R` | 14 | **0** |
| `pam.re.R` | 9 | **0** |
| `pam.rsph.metric.R` | 11 | **0** |
| `pam.rsph_metric.R` | 7 | **0** |
| `pam.schemper.R` | 10 | **0** |

Residual textual mentions in the filtered stream were investigated, not assumed
benign. Four mentions of `pam.rsph.R` and one of `r2nw` survive, all inside
**commit messages** such as *"Closes findings 38 and 39: R/pam.rsph.R was 92%
verbatim…"* and *"r2nw is commented 'unweighted measure' but is algebraically
identical to Re"*. These are the preserved discussion of §1, not code.

## 7. Post-rewrite verification

Two independent checks. **Neither may be a plain word search for `r_e` or
`rsph`** — that produces false positives on the preserved documentation and
commit messages and would fail a correct rewrite.

```bash
# 7.1 every recorded blob is unreachable
while read p b; do
  git cat-file -e "$b" 2>/dev/null && echo "STILL REACHABLE: $p $b"
done < ~/timemetric-purge-blobs-YYYYMMDD.txt
# expect: no output

# 7.2 no borrowed fragment survives in any blob under R/ or man/
git rev-list --all | while read c; do
  git grep -l -E 'tempi\.eventi|f\.assegna\.surv|newvare|incom\.sum|Re\.imp|r2nw|dRti' \
    "$c" -- R man 2>/dev/null
done
# expect: no output

# 7.3 the purge paths appear in no tree
git log --all --name-only --pretty=format: -- \
  R/pam.rsph.R R/pam.re.R R/pam.rsph.metric.R R/pam.rsph_metric.R \
  R/pam.schemper.R man/pam.rsph.Rd man/pam.re.Rd man/pam.rsph.metric.Rd \
  man/pam.rsph_metric.Rd man/pam.schemper.Rd | sort -u
# expect: no output

# 7.4 history is intact otherwise
git rev-list --all --count        # expect 106
git log --oneline -5              # messages preserved, hashes changed

# 7.5 the package still builds and tests green on the rewritten history
Rscript -e 'testthat::test_local()'   # expect 498 pass / 0 fail / 0 skip
R CMD build . && R CMD check --as-cran TimeMetric_0.2.0.tar.gz
```

**Why 7.2 is scoped to `R/` and `man/`:** tested against real history, an
unscoped scan false-positives on every candidate fragment but one —
`findings.md` and the audit quote the Italian identifiers as evidence, and the
test files name `dRti` and `my.survfit` as assertion targets. Scoped to the code
surface, all seven markers are clean. `my.survfit` and `n.incom` are excluded as
markers because they also occur in `R/pam.rsph_metric.R`, itself on the purge
list, so they cannot distinguish success from failure.

## 8. Commands that would update GitHub — **NOT EXECUTED**

> **None of the following has been run.** Each requires separate, explicit
> maintainer approval. They rewrite published history on a public repository.

```bash
# 8.1 re-add the remote that filter-repo removed
cd TimeMetric-rewrite.git
git remote add origin https://github.com/toz015/TimeMetric.git

# 8.2 inspect precisely what would change, before pushing anything
git push --mirror --dry-run origin

# 8.3 THE DESTRUCTIVE STEP - overwrites all published refs
git push --mirror --force origin

# 8.4 alternatively, branch by branch, with lease protection
git push --force-with-lease origin main
git push --force-with-lease origin joss-revision
```

Prerequisites before 8.3 or 8.4:

1. Both backups of §4 taken and **restore-tested** from the bundle.
2. All §7 checks pass on the rewritten history.
3. **The full-dependency CI job with `NOT_CRAN=true` confirmed green after the
   push and before CRAN submission** — see `removal-verification.md`. A plain
   `R CMD check` runs 441 of 498 assertions, so CI is what closes the gap.
4. Branch protection on `main` temporarily lifted, then restored.
5. Every collaborator warned in advance and instructed to re-clone.
6. Open PRs closed or recreated; they cannot survive the rewrite.
7. Forks deleted, or the blobs remain reachable through them.
8. GitHub Support asked to purge cached views of the removed blobs.

## 9. Deliberately not part of this plan

No tag, no GitHub release, and no CRAN submission. Those are separate decisions
taken after the rewrite has landed and CI is green.

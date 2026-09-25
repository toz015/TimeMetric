# History Purge Plan — combined borrowed-code and key-material purge

**Date prepared:** 2026-09-24
**Status:** **PREPARED ONLY. NOT EXECUTED.** No rewrite, push, force-push, tag,
release, or CRAN submission has been performed.
**Requires:** explicit maintainer approval before any command below is run.
**Verified with:** `git filter-repo --dry-run` in a throwaway mirror clone, which
writes export streams and rewrites nothing. `git-filter-repo` version
`a40bce548d2c`, at `/opt/homebrew/bin/git-filter-repo`. The bundle backup and
its restore were separately exercised end to end (§4).
**Revised 2026-09-24** on maintainer instruction: combined purge set defined
(§2), commit-count identity replaced as a success condition (§7.1), and
`push --mirror --force` demoted from the default publishing path (§8).

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

## 2. Final combined purge set

The borrowed-code purge and the key-material purge (Phase 1 spec step 8) are a
**single rewrite**. Two force-pushes would break every clone twice for no
benefit, so the sets are combined here and the history is rewritten once.

### 2.0 The combined set — 13 paths, 37 blobs

| Group | Path | Blobs |
|---|---|---|
| Key material | `q` | 1 |
| Key material | `q.pub` | 1 |
| Key material | `docs/security/key-exposure-report.md` | 1 |
| Borrowed code | `R/pam.rsph.R` | 11 |
| Borrowed code | `R/pam.re.R` | 4 |
| Borrowed code | `R/pam.rsph.metric.R` | 1 |
| Borrowed code | `R/pam.rsph_metric.R` | 4 |
| Borrowed code | `R/pam.schemper.R` | 6 |
| Borrowed docs | `man/pam.rsph.Rd` | 1 |
| Borrowed docs | `man/pam.re.Rd` | 3 |
| Borrowed docs | `man/pam.rsph.metric.Rd` | 1 |
| Borrowed docs | `man/pam.rsph_metric.Rd` | 1 |
| Borrowed docs | `man/pam.schemper.Rd` | 2 |

**On the key material:** the ed25519 pair was assessed and is **not a security
incident**. The maintainer confirmed it was never registered as a GitHub key, a
deploy key, or in any `authorized_keys`, so **no revocation is required**.
Removal is hygiene. Do not re-raise it as an incident.
`docs/security/key-exposure-report.md` is included because it is the report
about that material and was already moved out of the repository in `fc6c120`.

### 2.1 How the borrowed-code paths were derived

The ten borrowed paths were derived by scanning the **content** of every blob in
every commit for implementation-specific fragments, restricted to `R/` and
`man/`. The list was not copied from the design spec — and regenerating it mattered: **the spec listed 9
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
refs/heads/joss-revision    394f881   local only, NOT on the remote
refs/heads/main             1f54e6e   PUBLISHED - identical to origin/main
refs/remotes/origin/main    1f54e6e
refs/remotes/origin/HEAD    1f54e6e
tags                        (none)
```

* **107 commits** on all refs. Commits touching any path in the combined set,
  plus every descendant, get new hashes — in practice **all 107 are rewritten**.
* **Two commits are pruned entirely** (§7.1), leaving **105**.
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

# 4.3 Record the blobs that must become unreachable (COMBINED set)
for p in R/pam.rsph.R R/pam.re.R R/pam.rsph.metric.R R/pam.rsph_metric.R \
         R/pam.schemper.R man/pam.rsph.Rd man/pam.re.Rd \
         man/pam.rsph.metric.Rd man/pam.rsph_metric.Rd man/pam.schemper.Rd \
         q q.pub docs/security/key-exposure-report.md; do
  for c in $(git log --all --pretty='%H' -- "$p"); do
    b=$(git rev-parse --quiet --verify "$c:$p" 2>/dev/null) && echo "$p $b"
  done
done | sort -u > ~/timemetric-purge-blobs-$(date +%Y%m%d).txt
wc -l ~/timemetric-purge-blobs-*.txt    # expect 37 (34 borrowed + 3 key)
```

### Restoration procedure — **tested 2026-09-24**

Creating and restoring from the bundle was exercised end to end; the results
below are observed, not assumed.

```bash
git bundle create ~/tm.bundle --all
git bundle verify ~/tm.bundle        # => "The bundle records a complete history."

# RESTORE: use --mirror. A plain clone is NOT sufficient.
git clone --mirror ~/tm.bundle ~/tm-restored.git
```

**Tested result.** `git clone --mirror` restored all four refs and the full
history:

```
394f881 refs/heads/joss-revision
1f54e6e refs/heads/main
1f54e6e refs/remotes/origin/HEAD
1f54e6e refs/remotes/origin/main
commit count: 107
```

All three key-material blobs were confirmed still present in the restored
backup, which is what makes it a usable rollback:

```
present: 9607605caca745a1ef022168a1190ebfc6b9b819   (q)
present: 578c265dda4df22621d49beb92abc179b29b6ea7   (q.pub)
present: 0df72bcb136be4cd06bc0f2955f2c65b78325c8b   (key-exposure-report.md)
```

**Do not restore with a plain `git clone`.** That was tested too: it produced
only `refs/heads/joss-revision`, with `main` arriving as
`refs/remotes/origin/main` and no local `main` branch. The history is all there,
but the branch layout is wrong for a rollback push.

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
  --path man/pam.schemper.Rd \
  --path q \
  --path q.pub \
  --path docs/security/key-exposure-report.md
```

Add `--force` only if filter-repo refuses because the clone is not considered
fresh. Note that `filter-repo` **removes the `origin` remote** after rewriting,
by design, so the push in §8 re-adds it explicitly.

The key-material paths are already included above: this **is** the combined
rewrite, so history is rewritten once, not twice.

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

### 7.1 Commit pruning — the corrected success condition

**"The commit count remains 107" is wrong and must not be used.** The combined
purge makes two commits empty, and `git filter-repo` prunes empty commits by
default. Expecting an unchanged count would fail a correct rewrite.

The two commits were identified by testing every commit for whether its entire
change lies inside the combined purge set:

```bash
for c in $(git rev-list --all --no-merges); do
  paths=$(git show --pretty=format: --name-only "$c" | sed '/^$/d' | sort -u)
  [ -z "$paths" ] && continue
  outside=$(comm -23 <(echo "$paths") <(sort purge-set.txt))
  [ -z "$outside" ] && echo "$c"
done
```

**Exactly two commits are expected to disappear:**

| Commit | Subject | Entire content |
|---|---|---|
| `c873cb9` | original version | `q`, `q.pub` |
| `7eab69b` | docs: record the key exposure report | `docs/security/key-exposure-report.md` |

Both are path-only commits that add nothing else. **107 -> 105.**

Deliberately *not* pruned, and verified:

* `fc6c120` "docs: move the security report out of the repository" **survives**.
  It deletes `docs/security/key-exposure-report.md` but also modifies
  `docs/superpowers/pre-rewrite-review.md` and `tests/testthat/Rplots.pdf`, so
  only its deletion entry is dropped.
* The single merge commit `26fb014` touches no purge path.

**The success conditions are:**

1. **Only the two documented commits are pruned.**
2. **No unrelated commit is lost.**
3. **All surviving file contents outside the purge set are byte-identical.**

```bash
# 7.1a exactly the two expected commits were pruned, and nothing else.
# filter-repo writes .git/filter-repo/commit-map: "<old> <new>",
# with new = 40 zeros for a pruned commit.
awk '$2 ~ /^0+$/ {print $1}' .git/filter-repo/commit-map | sort > /tmp/pruned-actual.txt
printf '%s\n' \
  c873cb9... \
  7eab69b... | sort > /tmp/pruned-expected.txt   # use full 40-char SHAs
diff /tmp/pruned-expected.txt /tmp/pruned-actual.txt \
  && echo "PASS: only the documented commits were pruned" \
  || echo "FAIL: pruning differs from the plan -- STOP and investigate"

# 7.1b no unrelated commit lost: every surviving old commit maps to a new one
awk '$2 !~ /^0+$/' .git/filter-repo/commit-map | wc -l    # expect 105
git rev-list --all --count                                # expect 105

# 7.1c surviving content outside the purge set is byte-identical.
# Compare blob SHAs of all non-purge paths, old tree vs new tree.
BACKUP=~/tm-restored.git
PURGE_RE='(R/pam\.(rsph|re|rsph\.metric|rsph_metric|schemper)\.R|man/pam\.(rsph|re|rsph\.metric|rsph_metric|schemper)\.Rd|^q$|^q\.pub$|docs/security/key-exposure-report\.md)'
while read old new; do
  case "$new" in *[!0]*) : ;; *) continue ;; esac        # skip pruned
  a=$(git -C "$BACKUP" ls-tree -r --full-tree "$old" | grep -vE "$PURGE_RE" | shasum -a 256 | cut -d' ' -f1)
  b=$(git           ls-tree -r --full-tree "$new" | grep -vE "$PURGE_RE" | shasum -a 256 | cut -d' ' -f1)
  [ "$a" = "$b" ] || echo "CONTENT DIFFERS: $old -> $new"
done < .git/filter-repo/commit-map
# expect: no output. Comparing blob SHAs proves byte-identity, not just
# equal filenames.
```

### 7.2 Blob and fragment checks

Two further independent checks. **Neither may be a plain word search for `r_e` or
`rsph`** — that produces false positives on the preserved documentation and
commit messages and would fail a correct rewrite.

```bash
# 7.2a every recorded blob is unreachable (all 37)
while read p b; do
  git cat-file -e "$b" 2>/dev/null && echo "STILL REACHABLE: $p $b"
done < ~/timemetric-purge-blobs-YYYYMMDD.txt
# expect: no output

# 7.2b no borrowed fragment survives in any blob under R/ or man/
git rev-list --all | while read c; do
  git grep -l -E 'tempi\.eventi|f\.assegna\.surv|newvare|incom\.sum|Re\.imp|r2nw|dRti' \
    "$c" -- R man 2>/dev/null
done
# expect: no output

# 7.2c the purge paths appear in no tree
git log --all --name-only --pretty=format: -- \
  R/pam.rsph.R R/pam.re.R R/pam.rsph.metric.R R/pam.rsph_metric.R \
  R/pam.schemper.R man/pam.rsph.Rd man/pam.re.Rd man/pam.rsph.metric.Rd \
  man/pam.rsph_metric.Rd man/pam.schemper.Rd \
  q q.pub docs/security/key-exposure-report.md | sort -u
# expect: no output

# 7.4 see 7.1 below -- commit-count identity is NOT a valid success condition

# 7.2d the package still builds and tests green on the rewritten history
Rscript -e 'testthat::test_local()'   # expect 498 pass / 0 fail / 0 skip
R CMD build . && R CMD check --as-cran TimeMetric_0.2.0.tar.gz
```

**Why 7.2b is scoped to `R/` and `man/`:** tested against real history, an
unscoped scan false-positives on every candidate fragment but one —
`findings.md` and the audit quote the Italian identifiers as evidence, and the
test files name `dRti` and `my.survfit` as assertion targets. Scoped to the code
surface, all seven markers are clean. `my.survfit` and `n.incom` are excluded as
markers because they also occur in `R/pam.rsph_metric.R`, itself on the purge
list, so they cannot distinguish success from failure.

## 8. Publishing — **NOT EXECUTED**

> **Nothing below has been run.** Each step needs separate, explicit maintainer
> approval. These rewrite published history on a public repository.

**`git push --mirror --force` is deliberately NOT the default.** A mirror push
forces *every* ref in the local repo onto the remote and deletes any remote ref
absent locally — too blunt for a repository whose exact remote state must be
respected. Use explicit per-branch pushes with leases instead.

### 8.1 Inspect the live remote first — always

The remote is the authority on what is being overwritten, not local
remote-tracking refs, which can be stale.

```bash
git ls-remote --heads --tags origin
```

Observed 2026-09-24:

```
1f54e6ed8a8205975a942b50ebaa168e70ee0d27	refs/heads/main
```

**One branch, no tags.** `joss-revision` is **not** on the remote. If this
listing differs at execution time, **stop** — someone has pushed since, and
every lease value below is stale.

### 8.2 Re-add the remote

`filter-repo` removes `origin` by design after rewriting.

```bash
cd TimeMetric-rewrite.git
git remote add origin https://github.com/toz015/TimeMetric.git
git ls-remote --heads --tags origin     # confirm 8.1 still holds
```

### 8.3 `main` — rewriting already-published history

This is the destructive step. `main` exists on the remote at a known SHA, and
that SHA is what the lease is tied to, so the push **fails safely** if anyone
has pushed in the meantime.

```bash
# dry run first: shows exactly what would move
git push --dry-run --force-with-lease=refs/heads/main:1f54e6ed8a8205975a942b50ebaa168e70ee0d27 \
  origin main

# the real thing
git push --force-with-lease=refs/heads/main:1f54e6ed8a8205975a942b50ebaa168e70ee0d27 \
  origin main
```

Use `--force-with-lease=<ref>:<old-sha>` with the SHA written out, never bare
`--force-with-lease`: the bare form leases against the local remote-tracking ref,
which a `git fetch` silently updates, quietly destroying the protection.

### 8.4 `joss-revision` — a new remote branch, not a rewrite

`joss-revision` has never been pushed, so publishing it **creates** a branch and
overwrites nothing. It needs no force and no lease.

```bash
git push --dry-run origin joss-revision
git push origin joss-revision
```

**Open question for the maintainer:** should `joss-revision` be published at
all? It is the remediation branch. Options: publish it for review before
merging; merge it into `main` and publish only `main`; or keep it local until
the JOSS resubmission. Nothing here assumes an answer.

### 8.5 Tags — none exist

```
git ls-remote --tags origin     # => no output, 2026-09-24
git tag -l                      # => no output
```

There are no tags locally or remotely, so **no tag rewriting or force-pushing is
required**. If a tag is created later it must be made *after* the rewrite, from
the rewritten commits; a tag made beforehand would point at an obsolete SHA.

### 8.6 The mirror push, for reference only

```bash
# NOT the recommended path. Forces every ref and DELETES remote refs
# that are absent locally. Use only for a full restore from backup (section 4).
git push --mirror --force origin
```

Its legitimate use is rollback: pushing the mirror backup back over the remote.

### Prerequisites before 8.3 or 8.4

1. Both backups from §4 taken, and the **tested** restore procedure exercised.
2. All §7 checks pass, including the two-commit pruning condition in §7.1.
3. **The full-dependency CI job with `NOT_CRAN=true` confirmed green after the
   push and before CRAN submission** — see `removal-verification.md`. A plain
   `R CMD check` runs 441 of 498 assertions, so CI closes the gap.
4. `git ls-remote` re-checked immediately before the push; leases updated if the
   remote has moved.
5. Branch protection on `main` temporarily lifted, then restored.
6. Every collaborator warned in advance and told to re-clone.
7. Open PRs closed or recreated; they cannot survive the rewrite.
8. Forks deleted, or the blobs remain reachable through them.
9. GitHub Support asked to purge cached views of the removed blobs.

## 9. Deliberately not part of this plan

No tag, no GitHub release, and no CRAN submission. Those are separate decisions
taken after the rewrite has landed and CI is green.

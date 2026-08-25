# Public-release data hygiene

Do not commit real cluster paths, personal home directories, host addresses,
credentials, raw sequencing data or manuscript-restricted results.

Keep working configurations outside the repository. Commit only sanitized
examples that use placeholders. Before every public release, run:

```bash
python scripts/check_public_release.py --root .
```

The same check runs automatically in continuous integration. If a credential is
ever exposed, revoke it at the provider, remove it from Git history and create a
replacement before publishing the repository.

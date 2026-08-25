#!/usr/bin/env python3
"""Fail when a public-release tree contains private or sensitive material."""

from __future__ import annotations

import argparse
import ipaddress
import re
from pathlib import Path


SKIP_DIRECTORIES = {".git", ".venv", "__pycache__", ".pytest_cache"}
SELF_PATH = Path("scripts/check_public_release.py")
CODE_SUFFIXES = {".py", ".pyx", ".r", ".sh", ".bash", ".zsh", ".pl", ".pm"}
RAW_DATA_SUFFIXES = (
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".bam",
    ".cram",
    ".vcf",
    ".vcf.gz",
    ".bcf",
    ".fasta",
    ".fa",
    ".fna",
    ".cool",
    ".mcool",
    ".delta",
)
SENSITIVE_NAMES = re.compile(
    r"(?i)(?:^|[._-])(?:id_rsa|id_ed25519|credentials?|secrets?|private[_-]?key)(?:$|[._-])"
)
EMAIL = re.compile(r"\b[A-Z0-9._%+-]+@([A-Z0-9.-]+\.[A-Z]{2,})\b", re.IGNORECASE)
IPV4 = re.compile(r"(?<![0-9.])(?:\d{1,3}\.){3}\d{1,3}(?![0-9.])")
CJK = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff]")

TEXT_RULES = (
    (
        "private-storage-root",
        re.compile(r"(?i)(?<![A-Za-z0-9_])/(?:work\d*|scratch|gpfs|lustre)(?:/|$)"),
    ),
    (
        "personal-home-directory",
        re.compile(r"(?i)(?<![A-Za-z0-9_])/(?:home|Users)/[^/$<>{}\s]+/"),
    ),
    (
        "windows-user-directory",
        re.compile(r"(?i)\b[A-Z]:[\\/]Users[\\/][^\\/\s]+[\\/]"),
    ),
    (
        "private-key-material",
        re.compile(r"-----BEGIN (?:RSA |OPENSSH |EC |DSA )?PRIVATE KEY-----"),
    ),
    (
        "access-token-prefix",
        re.compile(r"\b(?:AKIA[0-9A-Z]{16}|gh[pousr]_[A-Za-z0-9_]{20,})\b"),
    ),
    (
        "credential-in-url",
        re.compile(r"(?i)(?:[?&](?:token|signature|sig|x-amz-signature)=|https?://[^/@\s]+:[^/@\s]+@)"),
    ),
    (
        "credential-assignment",
        re.compile(
            r"(?i)\b(?:password|passwd|api[_-]?key|access[_-]?token|client[_-]?secret)"
            r"\s*[:=]\s*[\"']?(?!\$|<|\{|REPLACE|YOUR_|EXAMPLE|NOT_SET)[A-Za-z0-9+/_.=-]{8,}"
        ),
    ),
    (
        "long-encoded-secret-candidate",
        re.compile(r"(?<![A-Za-z0-9+/])[A-Za-z0-9+/]{120,}={0,2}(?![A-Za-z0-9+/])"),
    ),
)


def decode(path: Path) -> str | None:
    data = path.read_bytes()
    if b"\x00" in data:
        return None
    for encoding in ("utf-8", "utf-8-sig"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            continue
    return None


def is_public_ip(candidate: str) -> bool:
    try:
        address = ipaddress.ip_address(candidate)
    except ValueError:
        return False
    return not (address.is_loopback or address.is_unspecified)


def iter_public_files(root: Path):
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(root)
        if any(part in SKIP_DIRECTORIES for part in rel.parts):
            continue
        yield path, rel


def scan(root: Path) -> list[str]:
    findings: list[str] = []
    for path, rel in iter_public_files(root):
        rel_posix = rel.as_posix()
        lower_name = path.name.lower()
        if SENSITIVE_NAMES.search(path.name) or lower_name.endswith((".pem", ".key", ".p12")):
            findings.append(f"{rel_posix}: sensitive-filename")
        if lower_name.endswith(RAW_DATA_SUFFIXES):
            findings.append(f"{rel_posix}: raw-data-file")
        if rel == SELF_PATH:
            continue
        text = decode(path)
        if text is None:
            continue
        for line_number, line in enumerate(text.splitlines(), start=1):
            for rule_name, pattern in TEXT_RULES:
                if pattern.search(line):
                    findings.append(f"{rel_posix}:{line_number}: {rule_name}")
            for match in EMAIL.finditer(line):
                domain = match.group(1).lower()
                if domain not in {"example.com", "example.org", "example.net", "example.invalid"}:
                    findings.append(f"{rel_posix}:{line_number}: email-address")
            for match in IPV4.finditer(line):
                if is_public_ip(match.group(0)):
                    findings.append(f"{rel_posix}:{line_number}: network-address")
            if path.suffix.lower() in CODE_SUFFIXES and CJK.search(line):
                findings.append(f"{rel_posix}:{line_number}: non-English-code-text")
    return findings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    root = args.root.resolve()
    findings = scan(root)
    if findings:
        print("Public-release audit failed:")
        for finding in findings:
            print(f"- {finding}")
        return 1
    file_count = sum(1 for _ in iter_public_files(root))
    print(f"Public-release audit passed for {file_count} files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

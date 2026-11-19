#!/usr/bin/env python3
"""Build a CRISPRTarget database from a FASTA file and write the data manager JSON.

A CRISPRTarget database consists of two directories, each containing a FASTA file
named ``<db-id>.fa`` together with a nucleotide BLAST database that was created with
``makeblastdb -parse_seqids -blastdb_version 4`` and a samtools index
(``<db-id>.fa.fai``):

- the target database, built from the sequences that should be searched, and
- the shuffled control database, which provides the background for the p-value
  calculation (built with ``esl-shuffle``).

Both directories are created with the CRISPRTarget helper scripts ``MakeUserDB.sh``
and ``MakeUserShuffledDB.sh``. They are written below ``--out-dir`` and are moved to
the Galaxy data manager directory by Galaxy according to ``data_manager_conf.xml``.
"""

import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import zipfile

# Public reference databases. The Zenodo record is version independent (concept DOI
# 10.5281/zenodo.14839268), the record id and the checksums pin the files that were
# curated by the CRISPRTarget authors (see https://github.com/davidchyou/CRISPRTarget).
ZENODO_RECORD = "21869551"
PUBLIC_DATABASES = {
    "vhdb": (
        "all_vhdb_full.fna.zip",
        "cd3e715adb9a748282ac9b095fb5d024",
        "viral and other sequences (small test database, ~49 MB download)",
    ),
    "phage": (
        "phage.fa.zip",
        "371b8a0ecba92f54afa4e2d0af4efe8e",
        "phage genomes (~660 MB download)",
    ),
    "plasmid": (
        "plasmid.fa.zip",
        "545a4b8763bfc1af6b97898bd6a01e3d",
        "plasmid sequences (~3.4 GB download)",
    ),
}

FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fas", ".fsa", ".ffn")
CHUNK = 1024 * 1024
PROGRESS = 100 * 1024 * 1024


def log(message):
    print(message, file=sys.stderr, flush=True)


def fail(message):
    raise SystemExit(f"ERROR: {message}")


def run(command):
    log("Running: " + " ".join(command))
    process = subprocess.run(command, capture_output=True, text=True)
    if process.stdout.strip():
        log(process.stdout.strip()[-2000:])
    if process.returncode != 0:
        log(process.stderr.strip()[-2000:])
        fail(f"command failed with exit code {process.returncode}: {' '.join(command)}")


def download(url, destination, expected_md5=None):
    """Download ``url`` to ``destination`` and check its md5 sum."""
    log(f"Downloading {url}")
    md5 = hashlib.md5()
    request = urllib.request.Request(
        url, headers={"User-Agent": "galaxy-data-manager-crisprtarget"}
    )
    with urllib.request.urlopen(request) as response, open(destination, "wb") as out:
        total = int(response.headers.get("Content-Length") or 0)
        done = 0
        reported = 0
        while True:
            chunk = response.read(CHUNK)
            if not chunk:
                break
            out.write(chunk)
            md5.update(chunk)
            done += len(chunk)
            if done - reported >= PROGRESS:
                reported = done
                if total:
                    log(f"  {done / 1e6:.0f} MB of {total / 1e6:.0f} MB")
                else:
                    log(f"  {done / 1e6:.0f} MB")
    if expected_md5 and md5.hexdigest() != expected_md5:
        fail(
            f"the md5 sum of {url} is {md5.hexdigest()}, expected {expected_md5}. "
            "The file was probably updated upstream - please report this to the tool developers."
        )
    log(f"Downloaded {done / 1e6:.0f} MB (md5 {md5.hexdigest()})")


def decompress(source, destination):
    with gzip.open(source, "rb") as src, open(destination, "wb") as dst:
        shutil.copyfileobj(src, dst)


def unpack_archive(archive, workdir, fasta_name):
    """Extract the largest FASTA file of a zip, gzip or plain file to ``workdir``."""
    with open(archive, "rb") as handle:
        magic = handle.read(4)
    output = os.path.join(workdir, fasta_name)
    if magic[:2] == b"PK":
        with zipfile.ZipFile(archive) as zf:
            members = [
                info
                for info in zf.infolist()
                if not info.is_dir()
                and not info.filename.startswith("__MACOSX")
                and info.filename.lower().endswith(FASTA_SUFFIXES + (".gz",))
            ]
            if not members:
                fail(f"the archive {os.path.basename(archive)} does not contain a FASTA file")
            member = max(members, key=lambda info: info.file_size)
            log(f"Extracting {member.filename} from {os.path.basename(archive)}")
            if member.filename.lower().endswith(".gz"):
                with zf.open(member) as src, gzip.open(output, "wb") as dst:
                    shutil.copyfileobj(src, dst)
            else:
                with zf.open(member) as src, open(output, "wb") as dst:
                    shutil.copyfileobj(src, dst)
        return output
    if magic[:2] == b"\x1f\x8b":
        log(f"Decompressing {os.path.basename(archive)}")
        decompress(archive, output)
        return output
    if os.path.basename(archive) != fasta_name:
        shutil.copyfile(archive, output)
        return output
    return archive


def fasta_size(path):
    """Return the number of bases in a FASTA file."""
    size = 0
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                size += len(line.strip())
    if not size:
        fail(f"{path} does not contain any sequence")
    return size


def build_database(fasta, out_dir, db_id, fasta_name):
    """Build the target and the shuffled control database, return their directories."""
    target_dir = os.path.join(out_dir, db_id)
    control_dir = os.path.join(out_dir, f"{db_id}_shuffled")

    run(["MakeUserDB.sh", fasta, target_dir])
    log("Shuffling the sequences for the control database (this can take hours for large databases)")
    run(["MakeUserShuffledDB.sh", fasta, control_dir])

    for directory in (target_dir, control_dir):
        for name in (fasta_name, f"{fasta_name}.fai", f"{fasta_name}.nin", f"{fasta_name}.nhr"):
            path = os.path.join(directory, name)
            if not os.path.isfile(path):
                listing = ", ".join(sorted(os.listdir(directory))) or "nothing"
                # MakeUserDB.sh discards the output of makeblastdb, so missing files
                # are the only indication for a failed database build.
                fail(
                    f"{path} was not created, the database could not be built (the directory "
                    f"contains: {listing})"
                )
    return target_dir, control_dir


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--public-database", choices=sorted(PUBLIC_DATABASES))
    source.add_argument("--url", help="URL of a FASTA file, optionally gzipped or zipped")
    source.add_argument("--fasta", help="path of a local FASTA file, optionally gzipped or zipped")
    parser.add_argument("--db-id", required=True, help="value (identifier) of the database")
    parser.add_argument("--name", required=True, help="name of the database")
    parser.add_argument("--dbsize", type=int, default=None, help="effective database size")
    parser.add_argument("--out-dir", required=True, help="directory for the built databases")
    parser.add_argument("--out-json", required=True, help="path of the data manager JSON file")
    args = parser.parse_args(argv)

    db_id = args.db_id.strip()
    if not db_id or not re.fullmatch(r"[A-Za-z0-9_.+-]+", db_id):
        fail(
            "The database identifier must not be empty and may only contain the characters "
            "A-Z, a-z, 0-9, _, ., + and -"
        )
    # The databases are named after the identifier, which also shows up in the sequence
    # identifiers that MakeUserDB.sh creates (SEQ_<number>|<fasta file name>).
    fasta_name = f"{db_id}.fa"

    workdir = tempfile.mkdtemp(prefix="crisprtarget_dm_")
    try:
        if args.public_database:
            file_name, md5, _ = PUBLIC_DATABASES[args.public_database]
            archive = os.path.join(workdir, file_name)
            download(
                f"https://zenodo.org/records/{ZENODO_RECORD}/files/{file_name}",
                archive,
                expected_md5=md5,
            )
            fasta = unpack_archive(archive, workdir, fasta_name)
        elif args.url:
            archive = os.path.join(workdir, "download")
            download(args.url, archive)
            fasta = unpack_archive(archive, workdir, fasta_name)
        else:
            # the file that Galaxy hands over has no meaningful extension, therefore the
            # compression is detected from the content (like for the files from URLs)
            fasta = unpack_archive(args.fasta, workdir, fasta_name)

        if os.path.basename(fasta) != fasta_name:
            renamed = os.path.join(workdir, fasta_name)
            shutil.move(fasta, renamed)
            fasta = renamed

        size = fasta_size(fasta)
        log(f"The database contains {size} bases")

        os.makedirs(args.out_dir, exist_ok=True)
        target_dir, control_dir = build_database(fasta, args.out_dir, db_id, fasta_name)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    dbsize = args.dbsize if args.dbsize else fasta_size(os.path.join(target_dir, fasta_name))
    entry = {
        "value": db_id,
        "name": args.name,
        "db_path": db_id,
        "ctrl_db": f"{db_id}_shuffled",
        "dbsize": str(dbsize),
    }
    with open(args.out_json, "w") as handle:
        json.dump({"data_tables": {"crispr_target_databases": [entry]}}, handle)
    log(f"Created the data manager JSON for the database {entry['value']} (dbsize {dbsize})")


if __name__ == "__main__":
    main()

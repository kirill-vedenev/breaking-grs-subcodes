"""Verify saved chunks and reconstruct a THK1 kernel and its derived panel."""
import argparse
from array import array
import hashlib
import json
from math import comb
from pathlib import Path
import shutil
import struct
import sys

from artifact import ROOT, TRIALS, sha
from validate_inputs import validate


def header(path):
    with Path(path).open("rb") as stream:
        raw = stream.read(20)
    if len(raw) != 20 or raw[:4] != b"THK1":
        raise ValueError("Expected a THK1 coefficient matrix")
    k, d, n, count = struct.unpack("<4I", raw[4:])
    if (k, d, n, count) != (30, 7, comb(30, 7), 256):
        raise ValueError("Unexpected kernel dimensions")
    if Path(path).stat().st_size != 20 + count * ((n + 7) // 8):
        raise ValueError("Kernel size does not match its header")
    return k, d, n, count


def verify_storage(trial):
    source = ROOT / "trials" / trial
    validate(source / "inputs")
    meta = json.loads((source / "kernel/storage.json").read_text())
    total = hashlib.sha256()
    size = 0
    first = b""
    for part in meta["chunks"]:
        path = source / "kernel" / part["file"]
        if path.parent != source / "kernel" or part["offset"] != size:
            raise ValueError("Invalid kernel chunk path or offset")
        if path.stat().st_size != part["bytes"] or sha(path) != part["sha256"]:
            raise ValueError(f"Kernel chunk mismatch: {path}")
        with path.open("rb") as stream:
            if not first:
                first = stream.read(20)
                stream.seek(0)
            for block in iter(lambda: stream.read(8 << 20), b""):
                total.update(block)
        size += part["bytes"]
    expected = b"THK1" + struct.pack("<4I", 30, 7, comb(30, 7), 256)
    if first != expected or size != 20 + 256 * ((comb(30, 7) + 7) // 8):
        raise ValueError("Invalid saved kernel header or size")
    if size != meta["bytes"] or total.hexdigest() != meta["sha256"]:
        raise ValueError("Concatenated kernel digest mismatch")
    if sha(source / "inputs/public.json") != meta["public_sha256"]:
        raise ValueError("Kernel coordinate basis mismatch")
    return meta


def panel_from_kernel(kernel, output, expected=None):
    import numpy as np
    k, d, n, count = header(kernel)
    output.mkdir()
    panel = output / "kernel-panel.u64le"
    words = 8
    packed = np.memmap(kernel, dtype=np.uint8, mode="r", offset=20,
                       shape=(count, (n + 7) // 8))
    with panel.open("wb") as stream:
        for start in range(0, n, 32768):
            end = min(n, start + 32768)
            bits = np.unpackbits(np.ascontiguousarray(packed[:, start // 8:(end + 7) // 8]),
                                 axis=1, bitorder="little")[:, :end - start].T
            raw = np.packbits(bits, axis=1, bitorder="little")
            padded = np.zeros((end - start, words * 8), dtype=np.uint8)
            padded[:, :raw.shape[1]] = raw
            stream.write(padded.tobytes())
    del packed
    monomials = output / "monomials-colex.u64"
    mask = (1 << d) - 1
    with monomials.open("wb") as stream:
        for start in range(0, n, 32768):
            values = array("Q")
            for _ in range(min(32768, n - start)):
                values.append(mask)
                low = mask & -mask
                above = mask + low
                mask = above | (((above ^ mask) >> 2) // low)
            if sys.byteorder != "little":
                values.byteswap()
            values.tofile(stream)
    if expected:
        for name, path in (("panel", panel), ("monomials", monomials)):
            record = expected[name]
            if path.stat().st_size != record["bytes"] or sha(path) != record["sha256"]:
                raise ValueError(f"Rebuilt {name} differs from the original")
    return panel, monomials


def materialize(trial, output, panel=False, kernel=None):
    source = ROOT / "trials" / trial
    validate(source / "inputs")
    meta = verify_storage(trial) if kernel is None else None
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"Use a fresh destination: {output}")
    if kernel is not None:
        kernel = Path(kernel).resolve()
        header(kernel)
    output.mkdir(parents=True)
    shutil.copytree(source / "inputs", output / "inputs")
    dest = output / "kernel.thk1"
    if kernel is not None:
        shutil.copyfile(kernel, dest)
    else:
        with dest.open("wb") as stream:
            for part in meta["chunks"]:
                with (source / "kernel" / part["file"]).open("rb") as chunk:
                    shutil.copyfileobj(chunk, stream, 8 << 20)
        if sha(dest) != meta["sha256"]:
            raise ValueError("Reconstructed kernel digest mismatch")
    header(dest)
    if panel:
        panel_from_kernel(dest, output / "polynomials",
                          meta["derived_files"] if meta else None)
        if meta:
            shutil.copyfile(source / "kernel/index.json", output / "polynomials/index.json")
        else:
            public = json.loads((source / "inputs/public.json").read_text())
            (output / "polynomials/index.json").write_text(json.dumps(
                dict(count=256, public_sha256=public["public_sha256"]), indent=2) + "\n")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trial", choices=TRIALS)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--panel", action="store_true")
    args = parser.parse_args()
    if args.panel and args.output is None:
        parser.error("--panel requires --output")
    if args.output:
        materialize(args.trial, args.output, args.panel)
        print(f"{args.trial}: reconstructed and verified")
    else:
        verify_storage(args.trial)
        print(f"{args.trial}: saved kernel verified (256 polynomials)")


if __name__ == "__main__":
    main()

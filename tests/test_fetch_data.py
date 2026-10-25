"""tools/fetch_data.py, without the network.

The fetch test streams a small local archive laid out like OpenMC's library
through a file:// URL.
"""
import hashlib
import io
import os
import sys
import tarfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))

import fetch_data  # noqa: E402


def _archive(path, files):
    with tarfile.open(path, "w:xz") as tar:
        for name, data in files.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))


def test_fetch_writes_only_the_requested_nuclides(tmp_path):
    archive = tmp_path / "lib.tar.xz"
    _archive(archive, {
        "endfb-viii.0-hdf5/photon/Fe.h5": b"photon",
        "endfb-viii.0-hdf5/neutron/Fe56.h5": b"iron",
        "endfb-viii.0-hdf5/neutron/Pb208.h5": b"lead",
        "endfb-viii.0-hdf5/neutron/U235.h5": b"uranium",
    })
    dest = tmp_path / "endfb" / "neutron"
    missing = fetch_data.fetch(["Fe56", "U235", "Xx1"], str(dest),
                               url=archive.as_uri())
    assert missing == ["Xx1"]
    assert sorted(os.listdir(dest)) == ["Fe56.h5", "U235.h5"]
    assert (dest / "U235.h5").read_bytes() == b"uranium"


def test_verify_flags_missing_and_different_files(tmp_path, capsys):
    (tmp_path / "Fe56.h5").write_bytes(b"iron")
    (tmp_path / "Pb208.h5").write_bytes(b"not lead")
    sums = {n: hashlib.sha256(d).hexdigest()
            for n, d in (("Fe56", b"iron"), ("Pb208", b"lead"), ("U235", b"u"))}
    assert fetch_data.verify(["Fe56", "Pb208", "U235"], str(tmp_path), sums) == ["Pb208", "U235"]
    out = capsys.readouterr().out
    assert "Fe56    ok" in out and "DIFFERS" in out and "missing" in out


def test_the_checksum_list_covers_every_set():
    sums = fetch_data.load_checksums()
    for group in fetch_data.SETS.values():
        assert set(group) <= set(sums)
    assert all(len(d) == 64 for d in sums.values())

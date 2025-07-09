#!/usr/bin/env python3
"""
Rucio File Checker Script
-------------------------

This script queries CMS Rucio for datasets, checks file availability through specific XRootD sites,
and writes accessible ROOT file paths to text files. It maintains a dynamic white/black list of
sites based on success or failure of file access.

======================================
🛠️ Required Python Packages
======================================

Make sure the following Python packages are installed in your user environment:

1. uproot
2. rucio-client
3. XRootD (Python bindings)
4. fsspec_xrootd (for filesystem-based access, optional but recommended)

To install them:

    python3 -m pip install --user uproot rucio-client xrootd fsspec fsspec_xrootd

To check if they are installed:

    python3 -c "import uproot; print('✔ uproot OK')"
    python3 -c "from rucio.client import Client; print('✔ rucio-client OK')"
    python3 -c "import XRootD.client; print('✔ XRootD OK')"
    python3 -c "import fsspec_xrootd; print('✔ fsspec_xrootd OK')"

If `fsspec_xrootd` fails to install via pip (common at CERN), use:

    pip install --user git+https://github.com/hsf/fsspec_xrootd.git

======================================
🔐 Grid Proxy (required for Rucio)
======================================

Before running, you **must** have a valid CMS VOMS proxy. Create it with:

    voms-proxy-init -voms cms -rfc --valid 168:0

To check the proxy is valid:

    voms-proxy-info -exists -valid 0:20 && echo "✔ Proxy valid" || echo "✘ Proxy invalid"

======================================
▶️ How to Run
======================================

Execute the script:

    python3 script_rucio_check.py

The script will:
  - Query the datasets defined in `dataset_map`
  - Validate XRootD access to files
  - Save valid file paths to `filelist_TTHH.txt`, `filelist_TTZZ.txt`, etc.

"""

import subprocess
import os
import uproot
import socket
from rucio.client import Client

def check_port(port):
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        sock.bind(("0.0.0.0", port))
        available = True
    except Exception:
        available = False
    sock.close()
    return available

def get_proxy_path() -> str:
    try:
        subprocess.run("voms-proxy-info -exists -valid 0:20", shell=True, check=True)
    except subprocess.CalledProcessError:
        raise Exception("VOMS proxy expired or non-existing. Run: `voms-proxy-init -voms cms -rfc --valid 168:0`")
    proxy = subprocess.check_output("voms-proxy-info -path", shell=True, text=True).strip()
    return proxy

# Rucio configuration path (CMS-specific)
if "RUCIO_HOME" not in os.environ:
    os.environ["RUCIO_HOME"] = "/cvmfs/cms.cern.ch/rucio/current"

def get_rucio_client(proxy=None) -> Client:
    if not proxy:
        proxy = get_proxy_path()
    return Client()

def write_to_file(filename, path):
    with open(filename, "a") as f:
        f.write(path + "\n")

def query_rucio(dataset, WL, BL, client, output_file, vetoT3=False):
    print(f"\n🧩 Querying Rucio for dataset: {dataset}")
    try:
        files = client.list_files(scope="cms", name=dataset)
    except Exception as e:
        print(f"❌ Failed to list files in dataset {dataset}: {e}")
        return

    total = 0
    ok = 0
    for filedata in files:
        total += 1
        filename = filedata['name']
        did = {"scope": "cms", "name": filename}
        try:
            replica_info = list(client.list_replicas([did]))[0]
        except Exception as e:
            print(f"⚠️ Could not retrieve replicas for {filename}: {e}")
            continue

        states = replica_info["states"]
        redirector = "root://xrootd-cms.infn.it/"
        chosen_site = None
        for site in states:
            if states[site] == "AVAILABLE" and "Tape" not in site:
                site_clean = site.replace("_Disk", "")
                if site_clean in WL:
                    redirector = f"root://xrootd-cms.infn.it//store/test/xrootd/{site_clean}/"
                    chosen_site = site_clean
                    break
                elif site_clean in BL or (vetoT3 and "T3" in site_clean):
                    continue
                else:
                    redirector = f"root://xrootd-cms.infn.it//store/test/xrootd/{site_clean}/"
                    chosen_site = site_clean

        full_path = redirector + filename
        print(f"📦 Testing file: {full_path}")
        try:
            with uproot.open(full_path) as file:
                file["Events"]
            write_to_file(output_file, full_path)
            WL.add(chosen_site)
            ok += 1
            print("   ✔ File OK")
        except Exception as e:
            BL.add(chosen_site)
            print(f"   ✘ File FAILED: {e}")

    print(f"\n✅ Dataset: {dataset}")
    print(f"   Total files checked: {total}")
    print(f"   Successfully accessed: {ok}")
    print(f"   Failed: {total - ok}")
    print(f"   Final White List: {WL}")
    print(f"   Final Black List: {BL}")

def main():
    client = get_rucio_client()
    WL = {"T2_DE_RWTH", "T2_IT_Rome", "T2_CH_CERN", "T2_ES_CIEMAT", "T1_IT_CNAF"}
    BL = {"T2_UK_London_IC"}

    dataset_map = {
        "filelist_TTHH.txt": "/TTHH-HHto4B_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM",
        "filelist_TTZZ.txt": "/TTZZ-ZZto4B_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM"
    }

    for output_file, dataset in dataset_map.items():
        if os.path.exists(output_file):
            os.remove(output_file)
        query_rucio(dataset, WL, BL, client, output_file)

if __name__ == "__main__":
    main()


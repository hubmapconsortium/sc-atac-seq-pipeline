#!/usr/bin/env python3
import json

genome_build_data = {"genome": "hg38", "annotations": {"source": "gencode", "version": 39}}

if __name__ == "__main__":
    with open("genome_build.json", "w") as f:
        json.dump(genome_build_data, f, indent=4)

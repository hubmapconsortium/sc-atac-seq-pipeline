#!/usr/bin/env python3
import json

genome_build_data = {"genome": "grch38", "annotations": {"source": "GENCODE", "version": 35}}

if __name__ == "__main__":
    with open("genome_build.json", "w") as f:
        json.dump(genome_build_data, file=f, indent=4)

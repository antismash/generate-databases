#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import (
    FeatureLocation,
    SeqFeature,
)
import json
import os
import sys
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Type,
)


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--dir', help="Directory containing the antiSMASH DB output folders")
    group.add_argument('--lof', help="List of files containing file names of cluster GenBank files to parse")

    args = parser.parse_args()

    file_list = []  # type: List[str]

    if args.dir is not None:
        for root, _, files in os.walk(args.dir):
            for filename in files:
                if filename.endswith(".gbk") and "final" in filename:
                    full_name = os.path.join(root, filename)
                    file_list.append(full_name)
                    print("found", full_name)
    else:
        with open(args.lof, 'r') as handle:
            file_list = list(map(str.strip, handle.readlines()))

    run(file_list)


class AsdbCluster:
    __slots__ = (
        'accession',
        'cluster_nr',
        'cluster_type',
        'contig_edge',
        'description',
        'loci',
        'minimal',
    )
    def __init__(self,
                 accession: str,
                 cluster_nr: int,
                 cluster_type: str,
                 contig_edge: bool,
                 description: str,
                 loci: List["AsdbLocus"],
                 minimal: bool = True,
                ) -> None:
        self.accession = accession
        self.cluster_nr = cluster_nr
        self.cluster_type = cluster_type
        self.contig_edge = contig_edge

        self.description = description
        # Kill useless, noisy " biosynthetic gene cluster" string
        if self.description.endswith(" biosynthetic gene cluster"):
            self.description = self.description[:-26]

        # TODO: Make sure locus tags are unique?
        self.loci = loci

        self.minimal = minimal

    @classmethod
    def from_biopython(cls: Type["AsdbCluster"],
                       record: SeqRecord,
                       idx: int,
                       minimal: bool):
        accession = record.annotations['accessions'][0]

        description = record.description

        cluster = record.features[idx]
        cluster_nr = None
        for note in cluster.qualifiers['note']:
            if note.startswith('Cluster number:'):
                cluster_nr = int(note[15:])

        cluster_type = cluster.qualifiers['product'][0]
        contig_edge = cluster.qualifiers['contig_edge'][0] == 'True'

        cds_features = []  # type: List[SeqFeature]
        for feature in record.features[idx-10:]:
            if feature.type != 'CDS':
                continue

            if cluster.location.start <= feature.location.start <= feature.location.end <= cluster.location.end:
                cds_features.append(feature)

            if feature.location.start > cluster.location.end:
                break

        if cluster_nr is None:
            raise RuntimeError("Invalid cluster entry in {}".format(accession))

        loci = []  # type: List[AsdbLocus]
        for feature in cds_features:
            loci.append(AsdbLocus.from_biopython(feature, accession, cluster_nr))

        return cls(accession, cluster_nr, cluster_type, contig_edge, description, loci, minimal)

    def write_fasta(self, handle):
        """Write all loci to handle in FASTA format."""
        for locus in self.loci:
            header = ">{l.accession}|{l.identifier}\n".format(l=locus)
            handle.write(header)
            handle.write(locus.sequence)
            handle.write('\n')

    def write_long_fasta(self, handle):
        """Write all loci to handle in FASTA format with long headers."""
        for locus in self.loci:
            header = ">{l.accession}|{l.cluster_nr}|{l.location}|{l.identifier}|{l.safe_annotation}\n".format(l=locus)
            handle.write(header)
            handle.write(locus.sequence)
            handle.write('\n')

    def write_long_fasta_legacy(self, handle):
        """Write all loci to handle in FASTA format with long headers in legacy format."""
        for locus in self.loci:
            strand = "+" if locus.location.strand == 1 else "-"
            start = locus.location.nofuzzy_start + 1
            header = (">{l.accession}|c{l.cluster_nr}|{start}-{l.location.nofuzzy_end}|{strand}|"
                      "{l.safe_locus_tag}|{l.safe_annotation}|{l.safe_accession}\n".format(l=locus, start=start, strand=strand))
            handle.write(header)
            handle.write(locus.sequence)
            handle.write('\n')

    def write_tabs_legacy(self, handle):
        """Write cluster info in ClusterBlast legacy tabular format."""
        elements = [
            self.accession,
            self.description,
            "c{}".format(self.cluster_nr),
            self.cluster_type,
            ";".join(map(lambda l: l.safe_locus_tag, self.loci)),
            ";".join(map(lambda l: l.identifier, self.loci)),
        ]
        handle.write("\t".join(elements) + "\n")

    def write_tabs(self, handle):
        """Write cluster info in ClusterBlast tabular format."""
        elements = [
            self.accession,
            self.description,
            str(self.cluster_nr),
            self.cluster_type,
            "minimal" if self.minimal else "full",
            "incomplete" if self.contig_edge else "complete",
            ";".join(map(lambda l: l.identifier, self.loci)),
        ]
        handle.write("\t".join(elements) + "\n")

    def to_dict(self):
        """Get a dictionary representation for JSON output."""
        self_dict = {
            "id": self.accession,
            "description": self.description,
            "cluster_nr": self.cluster_nr,
            "cluster_type": self.cluster_type,
            "analysis": "minimal" if self.minimal else "full",
            "truncated": self.contig_edge,
            "loci": list(map(lambda l: l.to_dict(), self.loci)),
        }
        return self_dict


class AsdbLocus:
    __slots__ = (
        'accession',
        'annotation',
        'cluster_nr',
        'gene_id',
        'location',
        'locus_tag',
        'protein_accession',
        'sequence',
    )

    def __init__(self,
                 location: FeatureLocation,
                 accession: str,
                 annotation: str,
                 cluster_nr: int,
                 sequence: str,
                 gene_id: Optional[str] = None,
                 locus_tag: Optional[str] = None,
                 protein_accession: Optional[str] = None) -> None:
        self.location = location
        self.accession = accession
        self.annotation = annotation
        self.cluster_nr = cluster_nr
        self.sequence = sequence
        self.gene_id = gene_id
        self.locus_tag = locus_tag
        self.protein_accession = protein_accession
        if not locus_tag and not gene_id and not protein_accession:
            raise RuntimeError("No valid identifier for feature")

        if not sequence:
            raise RuntimeError("No valid sequence")

    @classmethod
    def from_biopython(cls: Type["AsdbLocus"], feature: SeqFeature, accession:str, cluster_nr: int):
        annotation = feature.qualifiers.get('product', ["(unknown)"])[0]
        sequence = feature.qualifiers['translation'][0]

        optional_values = {}  # type: Dict[str, str]
        if 'gene' in feature.qualifiers:
            optional_values['gene_id'] = feature.qualifiers['gene'][0]
        if 'locus_tag' in feature.qualifiers:
            optional_values['locus_tag'] = feature.qualifiers['locus_tag'][0]
        if 'protein_id' in feature.qualifiers:
            optional_values['protein_accession'] = feature.qualifiers['protein_id'][0]
        return cls(feature.location, accession, annotation, cluster_nr, sequence, **optional_values)

    @property
    def identifier(self):
        if self.locus_tag:
            return self.locus_tag
        if self.protein_accession:
            return self.protein_accession
        if self.gene_id:
            return self.gene_id

    @property
    def safe_locus_tag(self):
        return self.locus_tag or "no_locus_tag"

    @property
    def safe_accession(self):
        return self.protein_accession or self.identifier

    @property
    def safe_annotation(self):
        ann = self.annotation.replace(" ", "_")
        if ann == "(unknown)":
            ann = ""
        return ann

    def to_dict(self):
        """Return a dict representation for JSON outputs."""
        self_dict = {
            "id": self.identifier,
            "annotation": self.annotation,
            "location": str(self.location),
        }
        if self.protein_accession:
            self_dict['protein_accession'] = self.protein_accession
        if self.locus_tag:
            self_dict['locus_tag'] = self.locus_tag
        if self.gene_id:
            self_dict['gene_id'] = self.gene_id

        return self_dict


def run(file_list: List[str]) -> None:
    clusters = []  # type: List[AsdbCluster]

    for gbk_file in file_list:
        log_file = os.path.join(os.path.dirname(gbk_file), 'log')

        if not os.path.exists(log_file):
            print("Missing log file for", gbk_file, file=sys.stderr)
            continue

        with open(log_file, "r") as handle:
            minimal = 'minimal' in handle.readline()

        records = SeqIO.parse(gbk_file, "genbank")
        for record in records:
            record.features.sort(key=lambda x: (x.location.nofuzzy_start, x.location.nofuzzy_end))
            for idx, feature in enumerate(record.features):
                if feature.type != "cluster":
                    continue
                cluster = AsdbCluster.from_biopython(record, idx, minimal)
                clusters.append(cluster)

    clusters.sort(key=lambda x: (x.accession, x.cluster_nr))

    with open("proteins.fasta", "w") as prots, \
         open("verbose_proteins.fasta", "w") as verbose_prots, \
         open("clusters.txt", "w") as tabs, \
         open("legacy_clusters.txt", "w") as legacy_tabs, \
         open("legacy_proteins.fasta", "w") as legacy_prots:
        for cluster in clusters:
            cluster.write_fasta(prots)
            cluster.write_long_fasta(verbose_prots)
            cluster.write_tabs(tabs)
            cluster.write_tabs_legacy(legacy_tabs)
            cluster.write_long_fasta_legacy(legacy_prots)

    with open("clusters.json", "w") as handle:
        json.dump(list(map(lambda x: x.to_dict(), clusters)), handle)


if __name__ == "__main__":
    main()

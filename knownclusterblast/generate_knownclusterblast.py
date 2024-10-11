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
from typing import (
    List,
    Optional,
    Type,
)

from mibig.converters.shared.mibig import MibigEntry
from mibig.converters.shared.mibig.common import CompletenessLevel
from mibig.converters.shared.mibig.biosynthesis.classes.base import SynthesisType
from mibig.converters.shared.mibig.biosynthesis.classes.ribosomal import Ribosomal

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', default=DATA_DIR, help="Directory containing the MIBiG antiSMASH run files")

    args = parser.parse_args()

    run(args)


class MibigCluster:
    __slots__ = (
        'cluster_type',
        'completeness',
        'description',
        'loci',
        'mibig_id',
    )
    def __init__(self,
                 mibig_id: str,
                 cluster_type: str,
                 completeteness: CompletenessLevel,
                 description: str,
                 loci: List["MibigLocus"],
                ) -> None:
        self.mibig_id = mibig_id
        self.cluster_type = cluster_type
        self.completeness = completeteness

        self.description = description
        # Kill useless, noisy " biosynthetic gene cluster" string
        if self.description.endswith(" biosynthetic gene cluster"):
            self.description = self.description[:-26]

        # TODO: Make sure locus tags are unique?
        self.loci = loci


    @classmethod
    def from_biopython(cls: Type["MibigCluster"],
                       mibig_id: str,
                       description: str,
                       record: SeqRecord,
                       cluster_type: str,
                       completeness: CompletenessLevel,
                      ):

        loci = []  # type: List[MibigLocus]
        for feature in record.features:
            if feature.type != "CDS":
                continue

            loci.append(MibigLocus.from_biopython(feature, mibig_id))

        return cls(mibig_id, cluster_type, completeness, description, loci)

    def write_fasta(self, handle):
        """Write all loci to handle in FASTA format."""
        for locus in self.loci:
            header = ">{l.accession}|{l.identifier}\n".format(l=locus)
            handle.write(header)
            handle.write(locus.sequence)
            handle.write('\n')

    def write_old_fasta(self, handle):
        """Write all loci to handle in old-style FASTA format."""
        for locus in self.loci:
            strand = "+" if locus.location.strand == 1 else "-"
            start = int(locus.location.start) + 1
            end = int(locus.location.end)
            header = (">{l.accession}|c1|{start}-{end}|{strand}|"
                      "{l.safe_locus_tag}|{l.safe_annotation}|{l.identifier}\n".format(l=locus, start=start, end=end, strand=strand))
            handle.write(header)
            handle.write(locus.sequence)
            handle.write('\n')

    def write_tabs(self, handle):
        """Write cluster info in ClusterBlast tabular format."""
        elements = [
            self.mibig_id,
            self.description,
            "c1",
            self.cluster_type,
            ";".join(map(lambda l: l.safe_locus_tag, self.loci)),
            ";".join(map(lambda l: l.identifier, self.loci)),
        ]
        handle.write("\t".join(elements) + "\n")

    def to_dict(self):
        """Get a dictionary representation for JSON output."""
        self_dict = {
            "id": self.mibig_id,
            "description": self.description,
            "cluster_type": self.cluster_type,
            "completeness": self.completeness.value,
            "loci": list(map(lambda l: l.to_dict(), self.loci)),
        }
        return self_dict


class MibigLocus:
    __slots__ = (
        'accession',
        'annotation',
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
                 sequence: str,
                 gene_id: Optional[str] = None,
                 locus_tag: Optional[str] = None,
                 protein_accession: Optional[str] = None) -> None:
        self.location = location
        self.accession = accession
        self.annotation = annotation
        self.sequence = sequence
        self.gene_id = gene_id
        self.locus_tag = locus_tag
        self.protein_accession = protein_accession
        if not locus_tag and not gene_id and not protein_accession:
            raise RuntimeError("No valid identifier for feature")

        if not sequence:
            raise RuntimeError("No valid sequence")

    @classmethod
    def from_biopython(cls: Type["MibigLocus"], feature: SeqFeature, accession: str):
        annotation = feature.qualifiers.get('product', ["(unknown)"])[0]
        sequence = feature.qualifiers['translation'][0]

        optional_values: dict[str, str] = {}
        if 'gene' in feature.qualifiers:
            optional_values['gene_id'] = feature.qualifiers['gene'][0]
        if 'locus_tag' in feature.qualifiers:
            optional_values['locus_tag'] = feature.qualifiers['locus_tag'][0]
        if 'protein_id' in feature.qualifiers:
            optional_values['protein_accession'] = feature.qualifiers['protein_id'][0]
        return cls(feature.location, accession, annotation, sequence, **optional_values)

    @property
    def identifier(self):
        if self.protein_accession:
            return self.protein_accession
        if self.locus_tag:
            return self.locus_tag
        if self.gene_id:
            return self.gene_id

    @property
    def safe_locus_tag(self):
        return self.locus_tag or self.identifier

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


def run(args):
    directories = []
    for name in os.listdir(args.data_dir):
        directory = os.path.join(args.data_dir, name)
        if os.path.isdir(directory):
            directories.append(directory)

    clusters = []  # type: List[MibigCluster]

    for directory in directories:
        mibig_id = os.path.basename(directory)
        json_file = os.path.join(directory, "annotations.json")
        if not os.path.exists(json_file):
            continue
        with open(json_file, 'r') as handle:
            raw_entry = open(json_file, 'r').read().strip()
        entry = MibigEntry.from_json(json.loads(raw_entry))
        cluster_type = "+".join(parse_bgc_types(entry))
        completeness = entry.completeness
        description = "/".join([c.name for c in entry.compounds])
        acc = entry.accession

        genbank_file = os.path.join(directory, "{}.gbk".format(acc))
        record = SeqIO.read(genbank_file, 'genbank')

        clusters.append(MibigCluster.from_biopython(mibig_id, description, record, cluster_type, completeness))

    clusters.sort(key=lambda x: x.mibig_id)

    with open("proteins.fasta", "w") as prots, open("legacy_proteins.fasta", "w") as legacy_prots, open("clusters.txt", "w") as tabs:
        for cluster in clusters:
            cluster.write_fasta(prots)
            cluster.write_old_fasta(legacy_prots)
            cluster.write_tabs(tabs)

    with open("clusters.json", "w") as handle:
        json.dump(list(map(lambda x: x.to_dict(), clusters)), handle)


def parse_bgc_types(entry: MibigEntry) -> list[str]:
    """Parse the BGC type."""
    types: list[str] = []
    for bgc_type in entry.biosynthesis.classes:
        full_type = bgc_type.class_name.value
        subclass = bgc_type.extra_info.subclass
        if subclass == "other":
            full_type += ":other"
            details = getattr(bgc_type.extra_info, "details", None)
            if details and not details.startswith("converted from"):
                full_type += f"({details})"
        elif bgc_type.class_name == SynthesisType.RIBOSOMAL and subclass == "RiPP":
            full_type += ":RiPP"
            assert isinstance(bgc_type.extra_info, Ribosomal)  # make mypy happy
            info = bgc_type.extra_info.ripp_type
            if info:
                full_type += f":{info}"
        elif subclass not in ("Unknown", None):
            full_type += f":{subclass}"
        types.append(full_type)

    return types


if __name__ == "__main__":
    main()

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

from mibig.converters.read.top import Everything

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', default=DATA_DIR, help="Directory containing the MIBiG antiSMASH run files")

    args = parser.parse_args()

    run(args)


class MibigCluster:
    __slots__ = (
        'cluster_type',
        'description',
        'loci',
        'mibig_id',
        'minimal',
    )
    def __init__(self,
                 mibig_id: str,
                 cluster_type: str,
                 description: str,
                 loci: List["MibigLocus"],
                 minimal: bool = True,
                ) -> None:
        self.mibig_id = mibig_id
        self.cluster_type = cluster_type

        self.description = description
        # Kill useless, noisy " biosynthetic gene cluster" string
        if self.description.endswith(" biosynthetic gene cluster"):
            self.description = self.description[:-26]

        # TODO: Make sure locus tags are unique?
        self.loci = loci

        self.minimal = minimal

    @classmethod
    def from_biopython(cls: Type["MibigCluster"],
                       mibig_id: str,
                       description: str,
                       record: SeqRecord,
                       cluster_type: str,
                       minimal: bool):

        loci = []  # type: List[MibigLocus]
        for feature in record.features:
            if feature.type != "CDS":
                continue

            loci.append(MibigLocus.from_biopython(feature, mibig_id))

        return cls(mibig_id, cluster_type, description, loci, minimal)

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
            start = locus.location.nofuzzy_start + 1
            header = (">{l.accession}|c1|{start}-{l.location.nofuzzy_end}|{strand}|"
                      "{l.safe_locus_tag}|{l.safe_annotation}|{l.identifier}\n".format(l=locus, start=start, strand=strand))
            handle.write(header)
            handle.write(locus.sequence)
            handle.write('\n')

    def write_tabs(self, handle):
        """Write cluster info in ClusterBlast tabular format."""
        elements = [
            self.mibig_id,
            self.description,
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
            "completeness": "partial" if self.minimal else "full",
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

        optional_values = {}  # type: Dict[str, str]
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

    print(directories)

    for directory in directories:
        mibig_id = os.path.basename(directory)
        json_file = os.path.join(directory, "{}.json".format(mibig_id))
        with open(json_file, 'r') as handle:
            raw_entry = open(json_file, 'r').read().strip()
        entry = Everything(json.loads(raw_entry))
        cluster_type = "+".join(parse_bgc_types(entry))
        minimal = entry.cluster.minimal
        description = "/".join([c.compound for c in entry.cluster.compounds])
        acc = entry.cluster.loci.accession

        genbank_file = os.path.join(directory, "{}.region001.gbk".format(acc))
        record = SeqIO.read(genbank_file, 'genbank')

        clusters.append(MibigCluster.from_biopython(mibig_id, description, record, cluster_type, minimal))

    clusters.sort(key=lambda x: x.mibig_id)

    with open("proteins.fasta", "w") as prots, open("legacy_proteins.fasta", "w") as legacy_prots, open("clusters.txt", "w") as tabs:
        for cluster in clusters:
            cluster.write_fasta(prots)
            cluster.write_old_fasta(legacy_prots)
            cluster.write_tabs(tabs)

    with open("clusters.json", "w") as handle:
        json.dump(list(map(lambda x: x.to_dict(), clusters)), handle)


def parse_bgc_types(entry):
    """Parse the BGC type."""
    types = entry.cluster.biosynthetic_class
    if "NRP" in types:
        subtype = parse_nrp_subtype(entry)
        if subtype:
            types.remove("NRP")
            types.append(subtype)
    if "Polyketide" in types:
        subtypes = parse_pks_subtypes(entry)
        if subtypes:
            types.remove("Polyketide")
            types.extend(subtypes)
    if "Saccharide" in types:
        subtype = parse_saccharide_subtype(entry)
        if subtype:
            types.remove("Saccharide")
            types.append(subtype)
    if "RiPP" in types:
        subtype = parse_ripp_subtype(entry)
        if subtype:
            types.remove("RiPP")
            types.append(subtype)
    if "Other" in types:
        subtype = parse_other_subtype(entry)
        if subtype:
            types.remove("Other")
            types.append(subtype)

    return types


def parse_nrp_subtype(entry):
    """Parse NRP subtype."""
    try:
        subtype = entry.cluster.nrp.subclass
    except AttributeError:
        return None
    if subtype in (None, "Unknown", "Other"):
        return None
    # I don't know any non-other lipopeptide in our dataset
    if subtype == "Other lipopeptide":
        subtype = "Lipopeptide"
    elif subtype == "Ca+-dependent lipopeptide":
        subtype = "Lipopeptide:Ca+-dependent lipopeptide"
    return "NRP:" + subtype


def parse_pks_subtypes(entry):
    """Parse PKS subtypes."""
    subtypes = []
    try:
        synthases = entry.cluster.polyketide.synthases
        for synthase in synthases:
            subtype = synthase.subclass or []
            subtypes.extend(subtype)
        return ["Polyketide:{} polyketide".format(t) for t in subtypes]
    except AttributeError:
        return []


def parse_saccharide_subtype(entry):
    """Parse Saccharide subtype."""
    try:
        subtype = entry.cluster.saccharide.subclass
    except AttributeError:
        return None
    if subtype is None:
        return None
    subtype = subtype.capitalize()
    if subtype in ("Unknown", "Other"):
        return None
    if subtype == "Hybrid/tailoring":
        subtype += " saccharide"
    return "Saccharide:" + subtype


def parse_ripp_subtype(entry):
    """Parse RiPPs subtype."""
    try:
        subtype = entry.cluster.ripp.subclass
    except AttributeError:
        return None
    if subtype in (None, "Unknown"):
        return None
    return "RiPP:" + subtype


def parse_other_subtype(entry):
    """Parse "Other" subtype."""
    try:
        subtype = entry.cluster.other.subclass
    except AttributeError:
        return None
    if subtype in (None, "Unknown"):
        return None
    return "Other:" + subtype


if __name__ == "__main__":
    main()

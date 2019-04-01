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

GENBANK_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "genbanks")
JSON_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "json")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genbank-dir', default=GENBANK_DIR, help="Directory containing the MIBiG GenBank files")
    parser.add_argument('--json-dir', default=JSON_DIR, help="Directory containing the MIBiG JSON files")

    args = parser.parse_args()

    run(args)


class MibigCluster:
    __slots__ = (
        'cluster_nr',
        'cluster_type',
        'description',
        'loci',
        'mibig_id',
        'minimal',
    )
    def __init__(self,
                 mibig_id: str,
                 cluster_nr: int,
                 cluster_type: str,
                 description: str,
                 loci: List["MibigLocus"],
                 minimal: bool = True,
                ) -> None:
        self.mibig_id = mibig_id
        self.cluster_nr = cluster_nr
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
                       record: SeqRecord,
                       cluster_type: str,
                       minimal: bool):
        mibig_id = record.annotations['accessions'][0]
        cluster_nr = record.annotations['sequence_version']

        description = record.description
        loci = []  # type: List[MibigLocus]
        for feature in record.features:
            if feature.type != "CDS":
                continue

            loci.append(MibigLocus.from_biopython(feature, mibig_id, cluster_nr))

        return cls(mibig_id, cluster_nr, cluster_type, description, loci, minimal)

    def write_fasta(self, handle):
        """Write all loci to handle in FASTA format."""
        for locus in self.loci:
            header = ">{l.accession}|{l.cluster_nr}|{l.identifier}\n".format(l=locus)
            handle.write(header)
            handle.write(locus.sequence)
            handle.write('\n')

    def write_old_fasta(self, handle):
        """Write all loci to handle in old-style FASTA format."""
        for locus in self.loci:
            strand = "+" if locus.location.strand == 1 else "-"
            start = locus.location.nofuzzy_start + 1
            header = (">{l.accession}|c{l.cluster_nr}|{start}-{l.location.nofuzzy_end}|{strand}|"
                      "{l.safe_locus_tag}|{l.safe_annotation}|{l.identifier}\n".format(l=locus, start=start, strand=strand))
            handle.write(header)
            handle.write(locus.sequence)
            handle.write('\n')

    def write_tabs(self, handle):
        """Write cluster info in ClusterBlast tabular format."""
        elements = [
            self.mibig_id,
            self.description,
            "c{}".format(self.cluster_nr),
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
            "cluster_nr": self.cluster_nr,
            "cluster_type": self.cluster_type,
            "completeness": "partial" if self.minimal else "full",
            "loci": list(map(lambda l: l.to_dict(), self.loci)),
        }
        return self_dict


class MibigLocus:
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
    def from_biopython(cls: Type["MibigLocus"], feature: SeqFeature, accession:str, cluster_nr: int):
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
        if self.protein_accession:
            return self.protein_accession
        if self.locus_tag:
            return self.locus_tag
        if self.gene_id:
            return self.gene_id

    @property
    def safe_locus_tag(self):
        return self.locus_tag or "no_locus_tag"

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
    genbank_files = [os.path.join(args.genbank_dir, filename) for filename in os.listdir(args.genbank_dir)]

    clusters = []  # type: List[MibigCluster]

    for gbk_file in genbank_files:
        accession, _ = os.path.splitext(os.path.basename(gbk_file))
        json_file = os.path.join(args.json_dir, accession + '.json')

        if not os.path.exists(json_file):
            print("Missing JSON file for", accession, file=sys.stderr)
            continue

        with open(json_file, "r") as handle:
            mibig_info = json.load(handle)

        cluster_type = get_cluster_type(mibig_info)
        minimal = 'minimal' in mibig_info['general_params']

        records = SeqIO.parse(gbk_file, "genbank")
        for record in records:
            cluster = MibigCluster.from_biopython(record, cluster_type, minimal)
            clusters.append(cluster)

    clusters.sort(key=lambda x: x.mibig_id)

    with open("proteins.fasta", "w") as prots, open("clusters.txt", "w") as tabs:
        for cluster in clusters:
            cluster.write_old_fasta(prots)
            cluster.write_tabs(tabs)

    with open("clusters.json", "w") as handle:
        json.dump(list(map(lambda x: x.to_dict(), clusters)), handle)


def get_cluster_type(mibig_info: Dict[str, Any]) -> str:
    cluster_type = ""
    biosyn_class = mibig_info['general_params']['biosyn_class']
    for cls in biosyn_class:
        cls = cls.lower()
        if cls == 'nrp':
            cls = 'nrps'
        elif cls == 'polyketide':
            cls = get_pks_subclass(mibig_info)
        elif cls == 'ripp':
            cls = get_ripp_subclass(mibig_info)

        cluster_type += cls + "-"

    return cluster_type[:-1]


def get_pks_subclass(mibig_info: Dict[str, Any]) -> str:
    if not 'Polyketide' in mibig_info['general_params']:
        print("No Polyketide section defined in", mibig_info['general_params']['mibig_accession'], file=sys.stderr)
        return "polyketide"

    if not 'pks_subclass' in mibig_info['general_params']['Polyketide']:
        print("No pks_subclass defined in", mibig_info['general_params']['mibig_accession'], file=sys.stderr)
        return "polyketide"

    subclasses = mibig_info['general_params']['Polyketide']['pks_subclass']
    new_cls = ""
    for scls in subclasses:
        scls = scls.lower()
        if scls.startswith("trans-at"):
            scls = "transatpks"
        elif scls.endswith("type i"):
            scls = "t1pks"
        elif scls.endswith("typei"):
            # yay, typo
            scls = "t1pks"
        elif scls.endswith("type ii"):
            scls = "t2pks"
        elif scls.endswith("type iii"):
            scls = "t3pks"
        elif scls.startswith("pufa"):
            scls = "pufa"
        elif scls == "other":
            scls = "polyketide"
        new_cls += scls + "+"

    if new_cls == "":
        print("Empty pks_subclass defined in", mibig_info['general_params']['mibig_accession'], file=sys.stderr)
        return "polyketide"

    return new_cls[:-1]


def get_ripp_subclass(mibig_info: Dict[str, Any]) -> str:
    if not 'RiPP' in mibig_info['general_params']:
        print("No RiPP section defined in", mibig_info['general_params']['mibig_accession'], file=sys.stderr)
        return "ripp"

    if not 'ripp_subclass' in mibig_info['general_params']['RiPP']:
        print("No ripp_subclass section defined in", mibig_info['general_params']['mibig_accession'], file=sys.stderr)
        return "ripp"

    subclass = mibig_info['general_params']['RiPP']['ripp_subclass'].lower()

    # fix misspellings
    if subclass == "lantipeptide":
        subclass = "lanthipeptide"
    elif subclass == "head-to-tailcyclized peptide":
        subclass = "head-to-tail cyclized peptide"
    elif subclass == "lassopeptide":
        subclass = "lasso peptide"
    elif subclass == "lap / microcin":
        subclass = "lap"

    return subclass


if __name__ == "__main__":
    main()

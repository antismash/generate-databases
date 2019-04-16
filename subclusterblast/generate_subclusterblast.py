#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
import json
from typing import (
    Dict,
    List,
    Type,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file', help="File containing the old style FASTA inputs")
    parser.add_argument('cluster_table', help="File containing the old style tabular cluster definition")

    args = parser.parse_args()

    run(args.fasta_file, args.cluster_table)


def run(fasta_file: str, cluster_table: str) -> None:
    proteins = {}  # type: Dict[str, Protein]
    clusters = []  # type: List[Cluster]

    raw_proteins = list(SeqIO.parse(fasta_file, "fasta", generic_protein))
    for raw_protein in raw_proteins:
        protein = Protein.from_biopython(raw_protein)
        assert protein.legacy_id not in proteins, "duplicate id: %s" % protein.legacy_id
        proteins[protein.legacy_id] = protein

    for line in open(cluster_table, "r"):
        cluster = Cluster.from_as4_clusterdesc(line, proteins)
        clusters.append(cluster)
    
    with open("legacy_clusters.txt", "w") as legacy_tabs, \
         open("legacy_proteins.fasta", "w") as legacy_prots:
        for cluster in clusters:
            cluster.write_tabs_legacy(legacy_tabs)
            cluster.write_fasta_legacy(legacy_prots)


class Protein:
    __slots__ = (
        'annotation',
        'cluster_accession',
        'cluster_number',
        'gene_name',
        'legacy_id',
        'location',
        'protein_id',
        'sequence',
        'unique_id',
    )

    def __init__(self,
                 annotation: str,
                 cluster_accession: str,
                 cluster_number: int,
                 gene_name: str,
                 location: FeatureLocation,
                 protein_id: str,
                 sequence: Seq,
                ) -> None:
        self.annotation = annotation
        self.cluster_accession = cluster_accession
        self.cluster_number = cluster_number
        self.gene_name = gene_name
        self.location = location
        self.protein_id = protein_id
        self.sequence = sequence
        self.unique_id = "{s.cluster_accession}_{s.cluster_number}_{s.gene_name}".format(s=self)
        self.legacy_id = "{s.cluster_accession}_{s.cluster_number}_{s.protein_id}".format(s=self)

    @classmethod
    def from_biopython(cls: Type["Protein"], rec: SeqRecord) -> Type["Protein"]:
        # use description because that seems to have the full FASTA header
        (
            cluster_acc_nr,
            cluster_nr_str,
            pos,
            strand,
            name,
            annotation,
            protein_id
        ) = rec.description.split('|')

        start, end = pos.split('-')
        strand = 1 if strand == '+' else -1
        location = FeatureLocation(int(start)-1, int(end), strand)

        annotation = annotation.replace("_", " ")

        cluster_accession, redundant_cluster_nr_str = cluster_acc_nr.rsplit("_", 1)
        cluster_number = int(cluster_nr_str.split("c")[-1])
        assert int(redundant_cluster_nr_str) == cluster_number, "cluster number mismatch: %s vs. %s" % (cluster_acc_nr, cluster_nr_str)

        return cls(annotation, cluster_accession, cluster_number, name, location, protein_id, rec.seq)

    @property
    def safe_annotation(self):
        return self.annotation.replace(" ", "_")


class Cluster:
    __slots__ = (
        'accession',
        'annotation',
        'cluster_number',
        'product',
        'proteins',
    )

    def __init__(self,
                 accession: str,
                 annotation: str,
                 cluster_number: int,
                 product: str,
                 proteins: List[Protein],
                ) -> None:
        self.accession = accession
        self.annotation = annotation
        self.cluster_number = cluster_number
        self.product = product
        self.proteins = proteins

    @classmethod
    def from_as4_clusterdesc(cls: Type["Cluster"], line: str, all_proteins: Dict[str, Protein]) -> Type["Cluster"]:
        (
            cluster_acc_nr,
            annotation,
            cluster_nr_str,
            product,
            protein_ids,
            protein_ids_repeat,
        ) = line.rstrip().split('\t')
        assert protein_ids == protein_ids_repeat, "protein_ids mismatch on %s" % cluster_acc_nr
        protein_ids = protein_ids.split(";")

        proteins = []  # type: List[Protein]
        for protein_id in protein_ids:
            if not protein_id:
                continue
            full_id = "{}_{}".format(cluster_acc_nr, protein_id)
            proteins.append(all_proteins[full_id])

        cluster_accession, redundant_cluster_nr_str = cluster_acc_nr.rsplit("_", 1)
        cluster_number = int(cluster_nr_str.split("c")[-1])
        assert int(redundant_cluster_nr_str) == cluster_number, "cluster number mismatch: %s vs. %s" % (cluster_acc_nr, cluster_nr_str)

        annotation = annotation.replace("_", " ")

        return cls(cluster_accession, annotation, cluster_number, product, proteins)

    def write_tabs_legacy(self, handle):
        """Write cluster info in SubClusterBlast legacy tabular format."""
        elements = [
            self.accession,
            self.annotation,
            "c{}".format(self.cluster_number),
            self.product,
            ";".join(map(lambda l: l.unique_id, self.proteins)),
            ";".join(map(lambda l: l.protein_id, self.proteins)),
        ]
        handle.write("\t".join(elements) + "\n")

    def write_fasta_legacy(self, handle):
        """Write proteins in SubClusterBlast legacy FASTA format."""
        for protein in self.proteins:
            start = protein.location.start + 1
            end = int(protein.location.end)
            strand = "+" if protein.location.strand == 1 else "-"
            header = (">{p.cluster_accession}|c{p.cluster_number}|{start}-{end}|"
                      "{strand}|{p.unique_id}|{p.safe_annotation}|"
                      "{p.protein_id}\n".format(p=protein, strand=strand, start=start, end=end))
            handle.write(header)
            handle.write(str(protein.sequence))
            handle.write("\n")


if __name__ == "__main__":
    main()

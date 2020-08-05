#!/usr/bin/env python3

import argparse
import json
import os
from typing import (
    Dict,
    Iterator,
    List,
    Optional,
    Type,
    Union,
)

import antismash
from antismash.common import secmet


DEFAULT_AS_OPTIONS = antismash.config.build_config(["--minimal"], modules=antismash.main.get_all_modules())


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--dir', help="Directory containing the antiSMASH DB output folders")
    group.add_argument('--lof', help="List containing file names of antiSMASH JSON output to parse")

    args = parser.parse_args()

    file_list = []  # type: List[str]

    if args.dir is not None:
        for root, _, files in os.walk(args.dir):
            for filename in files:
                if filename.endswith(".json"):
                    full_name = os.path.join(root, filename)
                    file_list.append(full_name)
                    print("found", full_name)
    else:
        with open(args.lof, 'r') as handle:
            file_list = list(map(str.strip, handle.readlines()))

    run(file_list)


class AsdbProtocluster:
    __slots__ = (
        'product',
        'start',
        'end',
    )
    def __init__(self, product, start, end):
        self.product = product
        self.start = start
        self.end = end

    def to_dict(self) -> Dict[str, Union[str, int]]:
        return {
            "product": self.product,
            "start": self.start,
            "end": self.end,
        }

    @classmethod
    def from_secmet(cls: Type["AsdbProtocluster"], protocluster: secmet.Protocluster) -> "AsdbProtocluster":
        return cls(protocluster.product, protocluster.location.start, protocluster.location.end)


class AsdbRegion:
    __slots__ = (
        'accession',
        'start',
        'end',
        'cluster_type',
        'contig_edge',
        'description',
        'loci',
        'minimal',
        'protoclusters',
    )

    def __init__(self,
                 accession: str,
                 start: int,
                 end: int,
                 products: str,
                 contig_edge: bool,
                 description: str,
                 protoclusters: List[AsdbProtocluster],
                 loci: List["AsdbLocus"],
                 minimal: bool = True,
                ) -> None:
        self.accession = accession
        self.cluster_type = products
        self.contig_edge = contig_edge
        self.start = start
        self.end = end
        self.description = description
        self.loci = loci
        self.minimal = minimal
        self.protoclusters = protoclusters

    @classmethod
    def from_secmet(cls: Type["AsdbRegion"],
                       record: secmet.Record,
                       region: secmet.Region,
                       minimal: bool = False):
        accession = record.annotations['accessions'][0]

        description = record.description

        if len(region.location) > 500000:
            print("Cluster too large:", record.id, region)
            return None

        loci = []  # type: List[AsdbLocus]
        for cds in region.cds_children:
            loci.append(AsdbLocus.from_secmet(cds, accession))

        protoclusters = []  # type: List[AsdbProtocluster]
        for protocluster in region.get_unique_protoclusters():
            protoclusters.append(AsdbProtocluster.from_secmet(protocluster))

        return cls(accession, region.location.start, region.location.end,
                   region.get_product_string(), region.contig_edge, description,
                   protoclusters, loci, minimal)

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
            header = ">{l.accession}|{l.location}|{l.identifier}|{l.safe_annotation}\n".format(l=locus)
            handle.write(header)
            handle.write(locus.sequence)
            handle.write('\n')

    def write_tabs(self, handle):
        """Write cluster info in ClusterBlast tabular format."""
        elements = [
            self.accession,
            self.description,
            "%s-%s" % (self.start, self.end),
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
            "start": self.start,
            "end": self.end,
            "cluster_type": self.cluster_type,
            "protoclusters": [pc.to_dict() for pc in self.protoclusters],
            "truncated": self.contig_edge,
            "loci": list(map(lambda l: l.to_dict(), self.loci)),
        }
        return self_dict


class AsdbLocus:
    __slots__ = (
        'accession',
        'annotation',
        'gene_id',
        'identifier',
        'location',
        'locus_tag',
        'protein_accession',
        'sequence',
    )

    def __init__(self,
                 location: secmet.locations.Location,
                 accession: str,
                 annotation: str,
                 sequence: str,
                 identifier: str,
                 gene_id: Optional[str] = None,
                 locus_tag: Optional[str] = None,
                 protein_accession: Optional[str] = None) -> None:
        self.location = location
        self.accession = accession
        self.annotation = annotation
        self.sequence = sequence
        self.gene_id = gene_id
        self.locus_tag = locus_tag
        self.identifier = identifier
        self.protein_accession = protein_accession
        if not locus_tag and not gene_id and not protein_accession:
            raise RuntimeError("No valid identifier for feature")

        if not sequence:
            raise RuntimeError("No valid sequence")

    @classmethod
    def from_secmet(cls: Type["AsdbLocus"], feature: secmet.CDSFeature, accession: str):
        annotation = feature.product or "(unknown)"

        optional_values = {}  # type: Dict[str, str]

        if feature.locus_tag:
            locus_tag = feature.locus_tag
            if locus_tag.find('allorf') > -1:
                print("making", locus_tag, "unique:", end=' ')
                locus_tag = "{}_{}".format(accession, locus_tag)
                print(locus_tag)
            optional_values['locus_tag'] = locus_tag

        optional_values['gene_id'] = feature.gene
        optional_values['protein_accession'] = feature.protein_id

        return cls(feature.location, accession, annotation, feature.translation,
                   feature.get_name(), **optional_values)

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


def regenerate(filename) -> Iterator[secmet.Record]:
    results = antismash.common.serialiser.AntismashResults.from_file(filename)
    for record, module_results in zip(results.records, results.results):
        record.strip_antismash_annotations()
        if antismash.detection.hmm_detection.__name__ not in module_results:
            return
        # reannotate in the order that antismash does, for correct regions etc
        antismash.main.run_detection(record, DEFAULT_AS_OPTIONS, module_results)

        def regen(raw, module):
            assert raw
            regenerated = module.regenerate_previous_results(raw, record, DEFAULT_AS_OPTIONS)
            assert regenerated is not None, "%s results failed to generate for %s" % (module.__name__, record.id)
            regenerated.add_to_record(record)
            return regenerated

        for module in antismash.main.get_analysis_modules():
            if module.__name__ in module_results:
                module_results[module.__name__] = regen(module_results[module.__name__], module)

        for val in module_results.values():
            assert not isinstance(val, dict)

        yield record


def run(file_list: List[str]) -> None:
    regions = []  # type: List[AsdbRegion]

    for filename in file_list:
        for record in regenerate(filename):
            for region in record.get_regions():
                regions.append(AsdbRegion.from_secmet(record, region))

    regions.sort(key=lambda x: (x.accession, x.start))

    with open("proteins.fasta", "w") as prots, \
         open("verbose_proteins.fasta", "w") as verbose_prots, \
         open("clusters.txt", "w") as tabs:
        for region in regions:
            region.write_fasta(prots)
            region.write_long_fasta(verbose_prots)
            region.write_tabs(tabs)

    with open("clusters.json", "w") as handle:
        json.dump(list(map(lambda x: x.to_dict(), regions)), handle)


if __name__ == "__main__":
    main()

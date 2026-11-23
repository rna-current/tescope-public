"""Known assembly registry.

TEscope downloads RepeatMasker annotation from UCSC goldenPath (the ``.fa.out.gz``
files are native RepeatMasker output, generated with the RepeatMasker/Repbase
runs UCSC documents per assembly). The original downloaded file is retained
next to the parsed cache for provenance. Redistribution note: TEscope never
ships these files; each user downloads them directly from the source.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class AssemblySpec:
    name: str
    organism: str
    description: str
    rmsk_url: str
    chromsizes_url: str
    source: str = "ucsc-goldenpath"

    def to_dict(self) -> dict:
        return asdict(self)


_UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath"

KNOWN_ASSEMBLIES: dict[str, AssemblySpec] = {
    "hg38": AssemblySpec(
        name="hg38",
        organism="Homo sapiens",
        description="Human GRCh38/hg38",
        rmsk_url=f"{_UCSC}/hg38/bigZips/hg38.fa.out.gz",
        chromsizes_url=f"{_UCSC}/hg38/bigZips/hg38.chrom.sizes",
    ),
    "hg19": AssemblySpec(
        name="hg19",
        organism="Homo sapiens",
        description="Human GRCh37/hg19",
        rmsk_url=f"{_UCSC}/hg19/bigZips/hg19.fa.out.gz",
        chromsizes_url=f"{_UCSC}/hg19/bigZips/hg19.chrom.sizes",
    ),
    "mm39": AssemblySpec(
        name="mm39",
        organism="Mus musculus",
        description="Mouse GRCm39/mm39",
        rmsk_url=f"{_UCSC}/mm39/bigZips/mm39.fa.out.gz",
        chromsizes_url=f"{_UCSC}/mm39/bigZips/mm39.chrom.sizes",
    ),
    "mm10": AssemblySpec(
        name="mm10",
        organism="Mus musculus",
        description="Mouse GRCm38/mm10",
        rmsk_url=f"{_UCSC}/mm10/bigZips/mm10.fa.out.gz",
        chromsizes_url=f"{_UCSC}/mm10/bigZips/mm10.chrom.sizes",
    ),
}

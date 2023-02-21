import re
import gzip
import sys
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord

def parse_gff3(file):
    """_reads in a gzipped gff3 file and returns a dictioary
    Input: 
        file: path to gzipped gff3 file
    Returns:
        dictionary: a dictionary of ENST:SeqRecords
    """

    # dicts
    ENST_info = dict()  # dict of seq records
    ENST_exon_dict = dict()  # a dict of dicts
    ENST_exon_dict2 = dict()  # a dict of req records
    ENST_CDS_dict = dict()
    ENST_codons_dict = dict()
    loc2exonID_dict = (
        dict()
    )  # the same location may have different ENSE IDs, and need ENST to distinguish
    loc2posType = (
        dict()
    )  # a dict that maps location to position types (e.g. 5UTR, exon intron, junction, 3UTR)


    # search patterns
    ENST_pattern = re.compile("(ENST.+?);")
    transcript_ID_pattern = re.compile("ID=transcript:(.+?);")
    Parent_pattern = re.compile("Parent=transcript:(.+?);")
    exon_id_pattern = re.compile("exon_id=(.+?);")
    CDS_id_pattern = re.compile("ID=CDS:(.+?);")
    rank_pattern = re.compile("rank=(.+?);")
    phase_pattern = re.compile("ensembl_phase=(.+?);")
    end_phase_pattern = re.compile("ensembl_end_phase=(.+?);")
    name_pattern = re.compile("Name=(.+?);")
    biotype_pattern = re.compile("biotype=(.+?);")
    version_pattern = re.compile("version=(.+?)")

    line_count = 0
    # go through gff3 file, store all ENST IDs
    print("-processing ENST IDs")
    with gzip.open(file, "rt") as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if len(fields) >= 2:
                chr = fields[0]
                type = fields[2]
                m = re.search(transcript_ID_pattern, line)  # require ID=ENSTXXXX
                if m:
                    ENST_id = m.group(1)
                    if ENST_id in ENST_info:
                        sys.exit(
                            f"{ENST_id} is already in ENST_info_dict"
                        )  # ID=transcript:ENSTXXXX should be unique
                    # get transcript parent gene and it's name
                    biotype = ""
                    version = "0"
                    name = ""
                    m2 = re.search(biotype_pattern, line)
                    m3 = re.search(version_pattern, line)
                    m4 = re.search(name_pattern, line)
                    if m2:
                        biotype = m2.group(1)
                    if m3:
                        version = m3.group(1)
                    if m4:
                        name = m4.group(1)
                    description = [f"{ENST_id}.{version}", biotype]
                    ENST_info[ENST_id] = SeqRecord(
                        "", id=ENST_id, description="|".join(description), name=name
                    )
                    line_count += 1

    # go through gff3 file, and store all exons in ENST_exon_dict
    print("-processing exons")
    with gzip.open(file, "rt") as fh:
        parent_ENST_id = ""
        for line in fh:
            # print(line.rstrip())
            fields = line.rstrip().split("\t")
            if len(fields) >= 2:
                chr = fields[0]
                type = fields[2]
                m2 = re.search(Parent_pattern, line)
                if m2 and type == "exon":
                    parent_ENST_id = m2.group(1)
                    if not parent_ENST_id in ENST_info.keys():
                        sys.exit(
                            f"{parent_ENST_id} is not found in ENST_info.keys(), offending line: {line.rstrip()} "
                        )
                    exon_id = (
                        exon_loc
                    ) = f"{parent_ENST_id}__{chr}_{fields[3]}-{fields[4]}_{fields[6]}"  # chr_start-end_strand
                    m3 = re.search(exon_id_pattern, line)
                    if m3:  # there is an exon ID
                        exon_id = m3.group(1)
                        loc2exonID_dict[exon_loc] = exon_id
                    exon_start = int(fields[3])
                    exon_end = int(fields[4])
                    exon_strand = plus_minus_strand_to_numeric(fields[6])
                    # append exon to seq record
                    ENST_info[parent_ENST_id].features.append(
                        SeqFeature(
                            location=FeatureLocation(
                                exon_start, exon_end, strand=exon_strand, ref=chr
                            ),
                            type=type,
                            id=exon_id,
                        )
                    )

    # parse GFF3 again, extracting CDS info, and referencing exon info
    print("-processing cds")
    with gzip.open(file, "rt") as fh:
        parent_ENST_id = ""
        for line in fh:
            # print(line.rstrip())
            fields = line.rstrip().split("\t")
            if len(fields) >= 2:
                chr = fields[0]
                type = fields[2]
                m2 = re.search(Parent_pattern, line)
                if m2 and type == "CDS":
                    parent_ENST_id = m2.group(1)
                    if not parent_ENST_id in ENST_info.keys():
                        sys.exit(
                            f"{parent_ENST_id} is not found in ENST_info.keys(), offending line: {line.rstrip()} "
                        )

                    CDS_id = f"{parent_ENST_id}__{chr}_{fields[3]}-{fields[4]}_{fields[6]}"  # chr_start-end_strand
                    # determine if CDS superimposes with an exon
                    if CDS_id in loc2exonID_dict.keys():
                        exon_id = loc2exonID_dict[CDS_id]
                        CDS_id = exon_id

                    # populate CDS info
                    CDS_start = int(fields[3])
                    CDS_end = int(fields[4])
                    CDS_strand = plus_minus_strand_to_numeric(fields[6])
                    CDS_phase = fields[7]
                    # append exon to seq record
                    # print(f"{CDS_start} {CDS_end} {CDS_strand} {type} {CDS_id}")
                    ENST_info[parent_ENST_id].features.append(
                        SeqFeature(
                            location=FeatureLocation(
                                CDS_start, CDS_end, strand=CDS_strand, ref=chr
                            ),
                            type=type,
                            id=CDS_id,
                        )
                    )
                    # copy over additional info from exon dict
                    # if parent_ENST_id in ENST_exon_dict.keys():
                    #     if exon_id in ENST_exon_dict[parent_ENST_id].keys():
                    #         ENST_CDS_dict[parent_ENST_id][CDS_id]["exon_rank"] = ENST_exon_dict[parent_ENST_id][exon_id]["rank"]
                    #         ENST_CDS_dict[parent_ENST_id][CDS_id]["exon_phase"] = ENST_exon_dict[parent_ENST_id][exon_id]["phase"]
                    #         ENST_CDS_dict[parent_ENST_id][CDS_id]["exon_end_phase"] = ENST_exon_dict[parent_ENST_id][exon_id]["end_phase"]

    # go through ENST_info, calculate span_start span_end for each ID
    for ID in ENST_info:
        if len(ENST_info[ID].features) > 0:
            coords = [int(i.location.end) for i in ENST_info[ID].features] + [
                int(i.location.start) for i in ENST_info[ID].features
            ]
            ENST_info[ID].span_start = min(coords)
            ENST_info[ID].span_end = max(coords)
            ENST_info[ID].chr = ENST_info[ID].features[0].ref

    return ENST_info

def plus_minus_strand_to_numeric(input):
    if input == "+":
        return 1
    elif input == "-":
        return -1
    else:
        return 0
from Bio.SeqFeature import FeatureLocation, ExactPosition

def extract_exon_intron_junc_coord(exons):
    """
    extract exon-intron junction coordinates
    input:  list of position-sorted exons, exon1,2,3... in the order that they are translated
    output: a list of exon-intron junction coordinates, 
            for each coordinate, the start position is the last base of the exon, the end position is the first base of the intron
    """
    assert(len(exons)>0)
    coord_list = []
    for idx, exon in enumerate(exons):
        if idx < len(exons) - 1: #skip last exon
            start = exon.location.start
            end = exon.location.end
            if exon.location.strand == 1:
                coord_list.append(FeatureLocation(ExactPosition(end), ExactPosition(end+1), strand=exon.location.strand, ref = exon.location.ref))
            else:
                coord_list.append(FeatureLocation(ExactPosition(start-1), ExactPosition(start), strand=exon.location.strand, ref = exon.location.ref))
    return(coord_list)

def extract_intron_exon_junc_coord(exons):
    """
    extract intron-exon junction coordinates
    input:  list of position-sorted exons, exon1,2,3... in the order that they are translated
    output: a list of intron-exon junction coordinates, 
            for each coordinate, the start position is the last base of the intron, the end position is the first base of the exon
    """
    assert(len(exons)>0)
    coord_list = []
    for idx, exon in enumerate(exons):
        if idx > 0: #skip first exon
            start = exon.location.start
            end = exon.location.end
            if exon.location.strand == 1:
                coord_list.append(FeatureLocation(ExactPosition(start-1), ExactPosition(start), strand=exon.location.strand, ref = exon.location.ref))
            else:
                coord_list.append(FeatureLocation(ExactPosition(end), ExactPosition(end+1), strand=exon.location.strand, ref = exon.location.ref))
    return(coord_list) 

def sort_exons(list_of_seqFeatures):
    if list_of_seqFeatures[0].strand == 1:
        sorted_list = sorted(list_of_seqFeatures, key=lambda seqFeature: seqFeature.location.start)
    else:
        sorted_list = sorted(list_of_seqFeatures, key=lambda seqFeature: seqFeature.location.start, reverse = True)
    return sorted_list


def extract_junc_seq(junc_coord, seq_dict, left_padding = 1, right_padding = 1):
    """
    extract junction sequence from a SeqRecord object
    input:  junction coordinates (FeatureLocation object)
            dictioary of seqRecords objects, indexed by record.id
    output: junction sequence (string)
    """
    seq_list = []
    for coord in junc_coord:

        seq_record = seq_dict[coord.ref]

        if coord.strand == 1:
            junc_seq = seq_record.seq[coord.start-left_padding:coord.end+right_padding-1]
            seq_list.append(junc_seq)
        else:
            junc_seq = seq_record.seq[coord.end-right_padding-1:coord.start+left_padding].reverse_complement()
            seq_list.append(junc_seq)
    return seq_list
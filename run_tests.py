from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition, SimpleLocation
import processGFF3 as pg
import processFa as pf
import utils as ut
from Bio.Seq import Seq

import unittest

class Tests_for_JuncSeq(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("loading genome sequence...")
        fasta_path = "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
        cls.genome_seqs = pf.parse_fasta(fasta_path)
        #print(cls.genome_seqs.keys())
    
        print("loading gff3 file...")
        gff3_path = "Homo_sapiens.GRCh38.107.gff3.gz"
        cls.ENST_info = pg.parse_gff3(gff3_path)


    def test_exon_intron(self):

        # -1 strand
        print("\ntesting extracting of exon-intron junctions with -1 stranded transcripts")
        features = self.ENST_info["ENST00000534980"].features
        exons = [i for i in features if i.type=="exon"]
        sorted_exons = ut.sort_exons(exons)
        junc_coord = ut.extract_exon_intron_junc_coord(sorted_exons)
        expected_result =  eval("[SimpleLocation(ExactPosition(118791097), ExactPosition(118791098), strand=-1, ref='11'),\
                             SimpleLocation(ExactPosition(118786051), ExactPosition(118786052), strand=-1, ref='11'),\
                             SimpleLocation(ExactPosition(118781120), ExactPosition(118781121), strand=-1, ref='11'),\
                             SimpleLocation(ExactPosition(118779631), ExactPosition(118779632), strand=-1, ref='11'),\
                             SimpleLocation(ExactPosition(118768222), ExactPosition(118768223), strand=-1, ref='11'),\
                             SimpleLocation(ExactPosition(118765208), ExactPosition(118765209), strand=-1, ref='11'),\
                             SimpleLocation(ExactPosition(118763211), ExactPosition(118763212), strand=-1, ref='11'),\
                             SimpleLocation(ExactPosition(118759921), ExactPosition(118759922), strand=-1, ref='11'),\
                             SimpleLocation(ExactPosition(118758773), ExactPosition(118758774), strand=-1, ref='11'),\
                             SimpleLocation(ExactPosition(118757170), ExactPosition(118757171), strand=-1, ref='11'),\
                             SimpleLocation(ExactPosition(118756259), ExactPosition(118756260), strand=-1, ref='11'),\
                             SimpleLocation(ExactPosition(118755401), ExactPosition(118755402), strand=-1, ref='11'),\
                             SimpleLocation(ExactPosition(118754704), ExactPosition(118754705), strand=-1, ref='11')]")
        self.assertEqual(junc_coord, expected_result)

        # +1 strand
        print("testing extracting of exon-intron junctions with +1 stranded transcripts")
        features = self.ENST_info["ENST00000380152"].features
        exons = [i for i in features if i.type=="exon"]
        sorted_exons = ut.sort_exons(exons)
        junc_coord = ut.extract_exon_intron_junc_coord(sorted_exons)
        expected_result =  eval("[SimpleLocation(ExactPosition(32315667), ExactPosition(32315668), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32316527), ExactPosition(32316528), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32319325), ExactPosition(32319326), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32325184), ExactPosition(32325185), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32326150), ExactPosition(32326151), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32326282), ExactPosition(32326283), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32326613), ExactPosition(32326614), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32329492), ExactPosition(32329493), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32331030), ExactPosition(32331031), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32333387), ExactPosition(32333388), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32341196), ExactPosition(32341197), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32344653), ExactPosition(32344654), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32346896), ExactPosition(32346897), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32355288), ExactPosition(32355289), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32356609), ExactPosition(32356610), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32357929), ExactPosition(32357930), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32362693), ExactPosition(32362694), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32363533), ExactPosition(32363534), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32370557), ExactPosition(32370558), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32371100), ExactPosition(32371101), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32376791), ExactPosition(32376792), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32379515), ExactPosition(32379516), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32379913), ExactPosition(32379914), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32380145), ExactPosition(32380146), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32394933), ExactPosition(32394934), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32397044), ExactPosition(32397045), strand=1, ref='13')]")
        self.assertEqual(junc_coord, expected_result)


    def test_intron_exon(self):

        # -1 strand
        print("\ntesting extracting of intron-exon junctions with -1 stranded transcripts")
        features = self.ENST_info["ENST00000534980"].features
        exons = [i for i in features if i.type=="exon"]
        sorted_exons = ut.sort_exons(exons)
        junc_coord = ut.extract_intron_exon_junc_coord(sorted_exons)
        expected_result =  eval("[SimpleLocation(ExactPosition(118786518), ExactPosition(118786519), strand=-1, ref='11'),\
                                    SimpleLocation(ExactPosition(118781184), ExactPosition(118781185), strand=-1, ref='11'),\
                                    SimpleLocation(ExactPosition(118779736), ExactPosition(118779737), strand=-1, ref='11'),\
                                    SimpleLocation(ExactPosition(118768352), ExactPosition(118768353), strand=-1, ref='11'),\
                                    SimpleLocation(ExactPosition(118765355), ExactPosition(118765356), strand=-1, ref='11'),\
                                    SimpleLocation(ExactPosition(118763306), ExactPosition(118763307), strand=-1, ref='11'),\
                                    SimpleLocation(ExactPosition(118760044), ExactPosition(118760045), strand=-1, ref='11'),\
                                    SimpleLocation(ExactPosition(118758902), ExactPosition(118758903), strand=-1, ref='11'),\
                                    SimpleLocation(ExactPosition(118757287), ExactPosition(118757288), strand=-1, ref='11'),\
                                    SimpleLocation(ExactPosition(118756323), ExactPosition(118756324), strand=-1, ref='11'),\
                                    SimpleLocation(ExactPosition(118755503), ExactPosition(118755504), strand=-1, ref='11'),\
                                    SimpleLocation(ExactPosition(118754887), ExactPosition(118754888), strand=-1, ref='11'),\
                                    SimpleLocation(ExactPosition(118752097), ExactPosition(118752098), strand=-1, ref='11')]")
        self.assertEqual(junc_coord, expected_result)

        # +1 strand
        print("testing extracting of intron-exon junctions with +1 stranded transcripts")
        features = self.ENST_info["ENST00000380152"].features
        exons = [i for i in features if i.type=="exon"]
        sorted_exons = ut.sort_exons(exons)
        junc_coord = ut.extract_intron_exon_junc_coord(sorted_exons)
        expected_result =  eval("[SimpleLocation(ExactPosition(32316421), ExactPosition(32316422), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32319076), ExactPosition(32319077), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32325075), ExactPosition(32325076), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32326100), ExactPosition(32326101), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32326241), ExactPosition(32326242), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32326498), ExactPosition(32326499), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32329442), ExactPosition(32329443), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32330918), ExactPosition(32330919), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32332271), ExactPosition(32332272), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32336264), ExactPosition(32336265), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32344557), ExactPosition(32344558), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32346826), ExactPosition(32346827), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32354860), ExactPosition(32354861), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32356427), ExactPosition(32356428), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32357741), ExactPosition(32357742), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32362522), ExactPosition(32362523), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32363178), ExactPosition(32363179), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32370401), ExactPosition(32370402), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32370955), ExactPosition(32370956), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32376669), ExactPosition(32376670), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32379316), ExactPosition(32379317), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32379749), ExactPosition(32379750), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32380006), ExactPosition(32380007), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32394688), ExactPosition(32394689), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32396897), ExactPosition(32396898), strand=1, ref='13'),\
                                    SimpleLocation(ExactPosition(32398161), ExactPosition(32398162), strand=1, ref='13')]")
        self.assertEqual(junc_coord, expected_result)


    def test_junc_seq_extract(self):
        
        # -1 strand
        print("\ntesting extracting of exon/intron junction sequences with -1 stranded transcripts")
        features = self.ENST_info["ENST00000534980"].features
        exons = [i for i in features if i.type=="exon"]
        sorted_exons = ut.sort_exons(exons)
        junc_coord = ut.extract_exon_intron_junc_coord(sorted_exons)
        res = ut.extract_junc_seq(junc_coord, self.genome_seqs, left_padding = 10, right_padding = 10)
        expected_res = eval("[Seq('AGCGGAGGAGGTAGGGCTGG'),\
                            Seq('CCACTATTAAGTAAGTGTTA'),\
                            Seq('CAAAACTTCGGTAAGTTGGT'),\
                            Seq('TCCTATTCAGGTATGTTAAC'),\
                            Seq('AATATACAAGGTTAGTTGAA'),\
                            Seq('GATGATACAGGTAAGAAATG'),\
                            Seq('ATTGGATGAGGTAATGTCTC'),\
                            Seq('GAAGTTCATGGTGAGTATAA'),\
                            Seq('TTTCTCCAGGGTAAGAAGTA'),\
                            Seq('AATGAGGCAGGTGAGTATAA'),\
                            Seq('GTTTGCACTGGTAAGTATTT'),\
                            Seq('GGAAGATCAGGTGAGGAAAA'),\
                            Seq('TAACAAGCATGTACGTCCCT')]")
        self.assertEqual(res, expected_res)

        print("\ntesting extracting of intron/exon junction sequences with -1 stranded transcripts")
        features = self.ENST_info["ENST00000534980"].features
        exons = [i for i in features if i.type=="exon"]
        sorted_exons = ut.sort_exons(exons)
        junc_coord = ut.extract_intron_exon_junc_coord(sorted_exons)
        res = ut.extract_junc_seq(junc_coord, self.genome_seqs, left_padding = 10, right_padding = 10)
        expected_res = eval("[Seq('TTTGTTTCAGATTGACGTGA'),\
                            Seq('TTCTCTCTAGACCTGGTGAT'),\
                            Seq('GTTTGACTAGGATGTGACCT'),\
                            Seq('TCTCTTTTAGGAGGAGAGCA'),\
                            Seq('ACTGTTTCAGCAATGGTGAT'),\
                            Seq('TCTTGCATAGTGCACGTGGT'),\
                            Seq('TCTACTGCAGGCAGATAAGT'),\
                            Seq('ccccccccAGAATTCCCATT'),\
                            Seq('CTGTGCTCAGCTTCAGATAA'),\
                            Seq('TTTCTCCTAGGAACATCGAA'),\
                            Seq('tttttttAAGATCTGTTTAC'),\
                            Seq('TGTTTTTTAGGTCGCTTTGG'),\
                            Seq('TCTTCTGTAGGCTTTGACAA')]")
        self.assertEqual(res, expected_res)

        # +1 strand
        print("\ntesting extracting of exon/intron junction sequences with +1 stranded transcripts")
        features = self.ENST_info["ENST00000380152"].features
        exons = [i for i in features if i.type=="exon"]
        sorted_exons = ut.sort_exons(exons)
        junc_coord = ut.extract_exon_intron_junc_coord(sorted_exons)
        res = ut.extract_junc_seq(junc_coord, self.genome_seqs, left_padding = 10, right_padding = 10)
        expected_res = eval("[Seq('TCTGGAGCGGGTTAGTGGTG'),\
                                Seq('AACAAAGCAGGTATTGACAA'),\
                                Seq('TTAGACTTAGGTAAGTAATg'),\
                                Seq('TTAGTGAAAGGTATGATGAA'),\
                                Seq('GATAAGTCAGGTATGATTAA'),\
                                Seq('GTTTGTGAAGGTAAATATTC'),\
                                Seq('GTGCTCATAGGTAATAATAG'),\
                                Seq('TACTACTGCTGTAAGTAAAT'),\
                                Seq('GCAAGTCATGGTAAGTCCTC'),\
                                Seq('GCTGATTCAGGTACCTCTGT'),\
                                Seq('ATCTTAGTGGGTAAGTGTTC'),\
                                Seq('ACTCCAGATGGTAAAATTAG'),\
                                Seq('TACCCTTTCGGTAAGACATG'),\
                                Seq('GAACCTTTAGGTATTGTATG'),\
                                Seq('TCATAAACAGGTATGTGTTT'),\
                                Seq('AATTTTATAGGTACTCTATG'),\
                                Seq('TAAAATACAGGCAAGTTTAA'),\
                                Seq('TATGTTAAAGGTAAATTAAT'),\
                                Seq('CCCTATACAGGTATGATGTA'),\
                                Seq('GAACATGAAGGTAAAATTAG'),\
                                Seq('TTACCTTGAGGTGAGAGAGT'),\
                                Seq('aaaGATTCAGGTAAGTATGT'),\
                                Seq('ACAACTACCGGTACAAACCT'),\
                                Seq('AAAAAAACAGGTAATGCACA'),\
                                Seq('TACTGTTGAGGTAAGGTTAC'),\
                                Seq('CAAGCTTCTGGTAAGTTAAT')]")
        self.assertEqual(res, expected_res)

        print("\ntesting extracting of intron/exon junction sequences with -1 stranded transcripts")
        features = self.ENST_info["ENST00000380152"].features
        exons = [i for i in features if i.type=="exon"]
        sorted_exons = ut.sort_exons(exons)
        junc_coord = ut.extract_intron_exon_junc_coord(sorted_exons)
        res = ut.extract_junc_seq(junc_coord, self.genome_seqs, left_padding = 10, right_padding = 10)
        expected_res = eval("[Seq('TGTTTTGCAGACTTATTTAC'),\
                            Seq('ttttAAATAGATTTAGGACC'),\
                            Seq('ACTGTTTCAGGAAGGAATGT'),\
                            Seq('TTTATTTTAGTCCTGTTGTT'),\
                            Seq('TTACCCCCAGTGGTATGTGG'),\
                            Seq('TTCCTCCCAGGGTCGTCAGA'),\
                            Seq('TATCTTACAGTCAGAAATGA'),\
                            Seq('ATTTTTGCAGAATGTGAAAA'),\
                            Seq('ACTTTAACAGGATTTGGAAA'),\
                            Seq('TTATGTTTAGGTTTATTGCA'),\
                            Seq('TTCTTTTTAGGAGAACCCTC'),\
                            Seq('TGTTTCCTAGGCACAATAAA'),\
                            Seq('CCCATTGCAGCACAACTAAG'),\
                            Seq('TCTTTGATAGATTTAATTAC'),\
                            Seq('ttttgtgtAGCTGTATACGT'),\
                            Seq('ATTTGTTCAGGGCTCTGTGT'),\
                            Seq('TCACTTTTAGATATGATACG'),\
                            Seq('ATTTGTCCAGATTTCTGCTA'),\
                            Seq('ATTATTACAGTGGATGGAGA'),\
                            Seq('GTTTTCTTAGAAAACACAAC'),\
                            Seq('ATGGTCACAGGGTTATTTCA'),\
                            Seq('CTCCAAACAGTTATACTGAG'),\
                            Seq('TTTTCTGTAGGTTTCAGATG'),\
                            Seq('tCCATTCTAGGACTTGCCCC'),\
                            Seq('ATTTTCTTAGAATATTGACA'),\
                            Seq('TTTTTATCAGATGTCTTCTC')]")
        self.assertEqual(res, expected_res)


 
if __name__ == '__main__':
    unittest.main()
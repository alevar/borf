import unittest

from classes.transcript import GTFObjectFactory, Object, Transcript, Exon, CDS
from utils.common import *

class TestObject(unittest.TestCase):
    def test_add_line_transcript(self):
        # Test adding a GTF line representing a transcript
        gtf_line = """chr1\tORFclust\ttranscript\t100\t200\t.\t+\t.\tgene_id "123"; transcript_id "456";"""
        obj = Object()
        result = obj.add_line(gtf_line)
        
        self.assertTrue(result)
        self.assertEqual(obj.get_type(), Types.Transcript)
        self.assertEqual(obj.get_seqid(), "chr1")
        self.assertEqual(obj.get_start(), 100)
        self.assertEqual(obj.get_end(), 201)  # Non-inclusive end
        self.assertEqual(obj.get_attr("gene_id"), "123")
        self.assertEqual(obj.get_attr("transcript_id"), "456")

    def test_add_line_exon(self):
        # Test adding a GTF line representing an exon
        gtf_line = """chr1\tORFclust\texon\t150\t180\t.\t+\t.\tgene_id "123"; transcript_id "456";"""
        obj = Object()
        result = obj.add_line(gtf_line)
        
        self.assertTrue(result)
        self.assertEqual(obj.get_type(), Types.Exon)
        self.assertEqual(obj.get_seqid(), "chr1")
        self.assertEqual(obj.get_start(), 150)
        self.assertEqual(obj.get_end(), 181)  # Non-inclusive end
        self.assertEqual(obj.get_attr("gene_id"), "123")
        self.assertEqual(obj.get_attr("transcript_id"), "456")

    def test_to_transcript(self):
        # Test converting the object to a Transcript object
        obj = Object()
        obj.set_seqid("chr1")
        obj.set_start(100)
        obj.set_end(200)
        obj.set_type(Types.Transcript)
        obj.add_attribute("gene_id", "123")
        obj.add_attribute("transcript_id", "456")
        
        transcript = obj.to_transcript()
        
        self.assertEqual(transcript.get_seqid(), "chr1")
        self.assertEqual(transcript.get_start(), 100)
        self.assertEqual(transcript.get_end(), 200)
        self.assertEqual(transcript.get_type(), Types.Transcript)
        self.assertEqual(transcript.get_attr("gene_id"), "123")
        self.assertEqual(transcript.get_attr("transcript_id"), "456")

    def test_overlaps(self):
        # Test checking if the object overlaps with another object
        obj1 = Object()
        obj1.set_start(100)
        obj1.set_end(200)
        
        obj2 = Object()
        obj2.set_start(150)
        obj2.set_end(250)
        
        obj3 = Object()
        obj3.set_start(300)
        obj3.set_end(400)
        
        self.assertTrue(obj1.overlaps(obj2))
        self.assertTrue(obj2.overlaps(obj1))
        self.assertFalse(obj1.overlaps(obj3))
        self.assertFalse(obj3.overlaps(obj1))

class TestTranscript(unittest.TestCase):
    def test_setup(self):
        gtf_line = """chr1\tORFclust\ttranscript\t100\t200\t.\t+\t.\tgene_id "123"; transcript_id "456";"""
        self.transcript = Transcript()
        self.transcript.add_line(gtf_line)
        self.exon1 = Exon()
        self.exon2 = Exon()
        self.exon3 = Exon()
        self.transcript.add_exon(self.exon1)
        self.transcript.add_exon(self.exon2)
        self.transcript.add_exon(self.exon3)

    def test_add_cds(self):
        gtf_line = """chr1\tORFclust\ttranscript\t100\t200\t.\t+\t.\tgene_id "123"; transcript_id "456";"""
        self.transcript = Transcript()
        self.transcript.add_line(gtf_line)
        cds = [(10, 20), (30, 40)]
        self.transcript.set_cds(cds)
        cds_regions = self.transcript.get_cds()
        print(cds_regions)
        self.assertEqual(len(cds_regions), 2)
        self.assertEqual(cds_regions[0][2].get_start(), 10)
        self.assertEqual(cds_regions[0][2].get_end(), 20)
        self.assertEqual(cds_regions[1][2].get_start(), 30)
        self.assertEqual(cds_regions[1][2].get_end(), 40)

if __name__ == '__main__':
    unittest.main()

import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from src.pairwise_alignment import perform_alignment

class TestPairwiseAlignment(unittest.TestCase):
    def setUp(self):
        self.seq1 = SeqRecord(Seq("ACGTACGT"), id="Seq1")
        self.seq2 = SeqRecord(Seq("ACGTTCGT"), id="Seq2")

    def test_global_alignment(self):
        alignment = perform_alignment(self.seq1, self.seq2, mode='global')
        self.assertIsNotNone(alignment)
        # Matches minus mismatch penalty
        self.assertTrue(alignment.score > 0) 

if __name__ == '__main__':
    unittest.main()
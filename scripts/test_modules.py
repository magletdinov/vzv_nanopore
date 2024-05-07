import unittest
from modules import *
import os
import numpy as np

unf8_map = {"A":65, "C":67, "G":71, "T": 84, "-": 45, "N": 78}

class TestAlignmenToNp(unittest.TestCase):
    def setUp(self):
        self.positive_control = "positive_control.fasta"
        with open(self.positive_control, 'w') as file:
            file.write(">ref\nACGTAAAA\n>seq1\nA-GTAAAT\n>seq2\nG-GTTAAA")
            
    def tearDown(self):
        os.remove(self.positive_control)
    
    def test_encode_and_align(self):
        expected_result = np.array([[unf8_map["A"], unf8_map["C"], unf8_map["G"], unf8_map["T"], unf8_map["A"], unf8_map["A"], unf8_map["A"], unf8_map["A"]],
                                    [unf8_map["A"], unf8_map["-"], unf8_map["G"], unf8_map["T"], unf8_map["A"], unf8_map["A"], unf8_map["A"], unf8_map["T"]],
                                    [unf8_map["G"], unf8_map["-"], unf8_map["G"], unf8_map["T"], unf8_map["T"], unf8_map["A"], unf8_map["A"], unf8_map["A"]]], dtype=np.int8)
        
        result = alignment_to_np(self.positive_control)
        
        self.assertTrue(np.array_equal(expected_result, result))
        self.assertEqual(expected_result.shape, result.shape)

        
class TestDeleteInsertions(unittest.TestCase):
    def test_delete_insertions(self):
        self.align_np = np.array([[65, 65, 45, 71, 67, 71, 45, 67],
                                  [65, 65, 84, 45, 71, 65, 65, 67],
                                  [45, 65, 84, 65, 67, 65, 65, 67]], dtype=np.int8)
        
        expected_result = np.array([[65, 65, 71, 67, 71, 67],
                                    [65, 65, 45, 71, 65, 67],
                                    [45, 65, 65, 67, 65, 67]], dtype=np.int8)
        
        result = delete_insertions(self.align_np)
        
        self.assertTrue(np.array_equal(expected_result, result))        
        
        
        
        
if __name__ == "__main__":
    unittest.main()      
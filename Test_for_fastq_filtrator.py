import unittest, os, sys
import random as rnd
from Bio.SeqUtils import GC
from Bio.SeqIO import parse

if not os.path.isfile("Fastq_filtrator.py"):
    print("Please, put \"Fastq_filtrator\" file in the directory with test-file!")
    exit()

from Fastq_filtrator import *

class fastQfiltratorTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open("my_fastq_file.fastq", "w") as test_fastQ_filter:
            letters = ("A", "G", "C", "T")
            test_file_str = ""
            quality = ""
            for i in range(1, 300):
                quality += "J"
                test_file_str += rnd.choice(letters)
                test_fastQ_filter.write(f"@line {i} \n")
                test_fastQ_filter.write(test_file_str + "\n")
                test_fastQ_filter.write("+" + "\n")
                test_fastQ_filter.write(quality + "\n")

    def test_file_name_parse(self):
        self.assertEqual(file_name_parse(("Fastq_filtrator.py", "my_fastq_file.fastq")),
                         "C:\\Users\\chera\\Desktop\\BI\\BI_2020_2021_Python\\BI_2020_2021_Python_fastq-filtrator-main\\my_fastq_file.fastq")
        self.assertEqual(file_name_parse(("Fastq_filtrator.py", "C:\\Users\\chera\\Desktop\\BI\\BI_2020_2021_Python\\BI_2020_2021_Python_fastq-filtrator-main\\my_fastq_file.fastq")),
                         "C:\\Users\\chera\\Desktop\\BI\\BI_2020_2021_Python\\BI_2020_2021_Python_fastq-filtrator-main\\my_fastq_file.fastq")

    def test_min_length_parse(self):
        self.assertEqual(min_length_parse(('--min_length', "20")), 20)
        self.assertRaises(SystemExit, min_length_parse, ('--min_length', "gan"))

    def test_outname_parse(self):
        self.assertEqual(outname_parse(("python", '--output_base_name', "filter_file", "my_fastq_file.fastq")), "filter_file")
        self.assertEqual(outname_parse(("python", "my_fastq_file.fastq")), "my_fastq_file")

    def test_keep_filtered_parse(self):
        self.assertEqual(keep_filtered_parse(("python",
                                              '--output_base_name', "filter_file",
                                              '--keep_filtered',
                                              "my_fastq_file.fastq")),
                         ["filter_file__passed.fastq", "filter_file__failed.fastq"])
        self.assertEqual(keep_filtered_parse(("python",
                                              '--output_base_name', "filter_file",
                                              "my_fastq_file.fastq")),
                         ["filter_file", "filter_file"])
        self.assertEqual(keep_filtered_parse(("python",
                                              '--keep_filtered',
                                              "my_fastq_file.fastq")),
                         ["my_fastq_file__passed.fastq", "my_fastq_file__failed.fastq"])
        self.assertEqual(keep_filtered_parse(("python",
                                              "my_fastq_file.fastq")),
                         ["my_fastq_file", "my_fastq_file"])

    def test_gc_boundry_parse(self):
        self.assertEqual(gc_boundry_parse(("python",
                                           '--gc_bounds',
                                           "30", "my_fastq_file.fastq")), [30])
        self.assertEqual(gc_boundry_parse(("python",
                                           '--gc_bounds', "30", "50",
                                           "my_fastq_file.fastq")), [30, 50])
        self.assertRaises(SystemExit, gc_boundry_parse, ("python",
                                           '--gc_bounds',
                                           "my_fastq_file.fastq"))

    def test_write_outfq(self):
        write_outfq("test_write_out_fq.txt", "test_write_outfq_function")
        with open("test_write_out_fq.txt", "r") as twf:
            first_string_file = twf.readline()
        self.assertEqual(first_string_file, "test_write_outfq_function\n")

    def test_filter_function(self):

        args1 = ("Fastq_filtrator.py", '--output_base_name', "filter_file_1", "my_fastq_file.fastq")
        filter_function(args1)
        self.assertEqual(os.path.isfile("filter_file_1.fastq"), True)


        args2 = ("Fastq_filtrator.py", '--min_length', "20", '--output_base_name', "filter_file_2", '--keep_filtered',
                 "my_fastq_file.fastq")
        filter_function(args2)
        self.assertEqual(os.path.isfile("filter_file_2__passed.fastq"), True)
        self.assertEqual(os.path.isfile("filter_file_2__failed.fastq"), True)
        with open("filter_file_2__passed.fastq") as fil_file:
            first_string = fil_file.readline()
        self.assertEqual(first_string.strip(), "@line 20")
        with open("filter_file_2__failed.fastq") as fil_file:
            first_string = fil_file.readline()
        self.assertEqual(first_string.strip(), "@line 1")


        args3 = ("Fastq_filtrator.py", '--min_length', "40", '--output_base_name', "filter_file_3", '--keep_filtered',
                 '--gc_bounds', "30", "60", "my_fastq_file.fastq")
        filter_function(args3)
        with open("filter_file_3__passed.fastq") as fil_file:
            first_string = fil_file.readline()
            second_string = fil_file.readline()
        self.assertEqual(len(second_string.strip()) >= 40, True)
        with open("filter_file_3__failed.fastq") as fil_file:
            first_string = fil_file.readline()
        self.assertEqual(first_string.strip(), "@line 1")
        for record in parse("filter_file_3__passed.fastq", "fastq-sanger"):
            self.assertEqual(GC(record.seq) >= 30, True)
            self.assertEqual(GC(record.seq) <= 60, True)
        for record in parse("filter_file_3__failed.fastq", "fastq-sanger"):
            if len(record.seq) > 39:
                self.assertEqual(GC(record.seq) > 30, False)
                self.assertEqual(GC(record.seq) < 60, False)


    @classmethod
    def tearDownClass(cls):
        os.remove("my_fastq_file.fastq")
        os.remove("test_write_out_fq.txt")
        os.remove("filter_file_1.fastq")
        os.remove("filter_file_2__passed.fastq")
        os.remove("filter_file_2__failed.fastq")
        os.remove("filter_file_3__passed.fastq")
        os.remove("filter_file_3__failed.fastq")

if __name__ == '__main__':
    unittest.main()
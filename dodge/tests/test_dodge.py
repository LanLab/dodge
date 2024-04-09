import unittest

from dodge.dodge import ident_meta_columns


class TestMetadatacols(unittest.TestCase):
    def test_correct_metacols(self):
        infile = open("dodge/tests/inputs/correct_mgt_metadata.txt", "r")
        inf = infile.read().splitlines()[0]
        infile.close()
        header = inf.split("\t")
        ident_meta_columns_out = ident_meta_columns(header)
        anymissing = any([x == "" for x in ident_meta_columns_out[:4]])
        self.assertFalse(anymissing, 'Isolate, Date, Month or Year column detection not working')

    def test_noisolate_metacols(self):
        infile = open("dodge/tests/inputs/noisolate_mgt_metadata.txt", "r")
        inf = infile.read().splitlines()[0]
        infile.close()
        header = inf.split("\t")
        # ident_meta_columns_out = ident_meta_columns(header)
        # anymissing = any([x == "" for x in ident_meta_columns_out[:4]])
        # self.assertFalse(anymissing, 'Isolate, Date, Month or Year column detection not working')
        self.assertRaises(SystemExit,ident_meta_columns,header,isolatename=False,enterobase=False)

    def test_empty_noisolate_metacols(self):
        infile = open("dodge/tests/inputs/empty_mgt_metadata.txt", "r")
        inf = infile.read().splitlines()[0]
        infile.close()
        header = inf.split("\t")
        self.assertRaises(SystemExit,ident_meta_columns,header,isolatename=False,enterobase=False)


if __name__ == '__main__':
    unittest.main()
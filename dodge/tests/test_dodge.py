import unittest
from datetime import datetime
from dodge.dodge import ident_meta_columns,get_start_end

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

class TestGetStartEnd(unittest.TestCase):

    def test_week_with_start_end_date(self):
        startdate, enddate = get_start_end("week", "2024-01-01", "2024-01-07", {}, {})
        self.assertEqual(startdate, datetime(2024, 1, 1))
        self.assertEqual(enddate, datetime(2024, 1, 7))

    def test_week_with_start_date_only(self):
        startdate, enddate = get_start_end("week", "2024-01-01", None, {}, {})
        self.assertEqual(startdate, datetime(2024, 1, 1))
        # End date should be the newest strain date in this case

    def test_week_with_end_date_only(self):
        startdate, enddate = get_start_end("week", None, "2024-01-07", {}, {})
        # Start date should be the oldest strain date in this case
        self.assertEqual(enddate, datetime(2024, 1, 7))

    def test_week_with_no_dates(self):
        startdate, enddate = get_start_end("week", None, None, {}, {})
        # Start date should be the oldest strain date and end date should be the newest strain date
        # But it depends on how `get_oldest_strain_date` and `get_newest_strain_date` are implemented

    def test_month_with_start_end_date(self):
        startdate, enddate = get_start_end("month", "2024-01", "2024-02", {}, {})
        self.assertEqual(startdate, datetime(2024, 1, 1))
        self.assertEqual(enddate, datetime(2024, 2, 1))


if __name__ == '__main__':
    unittest.main()
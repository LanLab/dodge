import os
import shutil
import unittest
from datetime import datetime
from dodge import dodge
import argparse
import hashlib
from time import sleep as sl


def compare_files(file1_path, file2_path):
    # Initialize SHA-256 hash objects
    hash1 = hashlib.sha256()
    hash2 = hashlib.sha256()

    # Open the first file and update hash object with its contents
    with open(file1_path, "rb") as file1:
        while True:
            chunk = file1.read(4096)
            if not chunk:
                break
            hash1.update(chunk)

    # Open the second file and update hash object with its contents
    with open(file2_path, "rb") as file2:
        while True:
            chunk = file2.read(4096)
            if not chunk:
                break
            hash2.update(chunk)

    # Compare the hashes
    return hash1.digest() == hash2.digest()


class TestMetadatacols(unittest.TestCase):
    def test_correct_metacols(self):
        cwd = os.getcwd()
        infile = open(cwd+"/dodge/tests/inputs/correct_mgt_metadata.txt", "r")
        inf = infile.read().splitlines()[0]
        infile.close()
        header = inf.split("\t")
        ident_meta_columns_out = dodge.ident_meta_columns(header)
        anymissing = any([x == "" for x in ident_meta_columns_out[:4]])
        self.assertFalse(anymissing, 'Isolate, Date, Month or Year column detection not working')

    def test_noisolate_metacols(self):
        cwd = os.getcwd()
        infile = open(cwd + "/dodge/tests/inputs/noisolate_mgt_metadata.txt", "r")
        inf = infile.read().splitlines()[0]
        infile.close()
        header = inf.split("\t")
        self.assertRaises(SystemExit,dodge.ident_meta_columns,header,isolatename=False,enterobase=False)

    def test_empty_noisolate_metacols(self):
        cwd = os.getcwd()
        infile = open(cwd + "/dodge/tests/inputs/empty_mgt_metadata.txt", "r")
        inf = infile.read().splitlines()[0]
        infile.close()
        header = inf.split("\t")
        self.assertRaises(SystemExit,dodge.ident_meta_columns,header,isolatename=False,enterobase=False)

class TestGetStartEnd(unittest.TestCase):

    def test_week_with_start_end_date(self):
        startdate, enddate = dodge.get_start_end("week", "2024-01-01", "2024-01-07", {}, {})
        self.assertEqual(startdate, datetime(2024, 1, 1))
        self.assertEqual(enddate, datetime(2024, 1, 7))

    def test_week_with_start_date_only(self):
        startdate, enddate = dodge.get_start_end("week", "2024-01-01", None, {}, {})
        self.assertEqual(startdate, datetime(2024, 1, 1))
        # End date should be the newest strain date in this case

    def test_week_with_end_date_only(self):
        startdate, enddate = dodge.get_start_end("week", None, "2024-01-07", {}, {})
        # Start date should be the oldest strain date in this case
        self.assertEqual(enddate, datetime(2024, 1, 7))


    def test_month_with_start_end_date(self):
        startdate, enddate = dodge.get_start_end("month", "2024-01", "2024-02", {}, {})
        self.assertEqual(startdate, datetime(2024, 1, 1))
        self.assertEqual(enddate, datetime(2024, 2, 1))

class TestSNPDistMetric(unittest.TestCase):

    def test_snp_dist_metric_no_missmatch(self):
        a = "ACGTGTAC"
        b = "ACGTGTAC"
        args = argparse.Namespace(max_missmatch=3)
        result = dodge.snp_dist_metric(a, b, args)
        self.assertEqual(result, 0)

    def test_snp_dist_metric_with_missmatch_within_threshold(self):
        a = "ACGTCTAC"
        b = "ACGTGTAC"
        args = argparse.Namespace(max_missmatch=3)
        result = dodge.snp_dist_metric(a, b, args)
        self.assertEqual(result, 1)

    def test_snp_dist_metric_with_missmatch_exceeding_threshold(self):
        a = "ACGATACC"
        b = "ACGTGTAC"
        args = argparse.Namespace(max_missmatch=3)
        result = dodge.snp_dist_metric(a, b, args)
        self.assertEqual(result, 3)  # Missmatch exceeds the threshold

    def test_snp_dist_metric_with_missing_values(self):
        a = "ACG-TACC"
        b = "ACGTGNAC"
        args = argparse.Namespace(max_missmatch=3)
        result = dodge.snp_dist_metric(a, b, args)
        self.assertEqual(result, 2)

class TestDodge_full(unittest.TestCase):
    def setUp(self):
        self.outputfolder = "outputs"
        if os.path.exists(self.outputfolder):
            shutil.rmtree(self.outputfolder)
            os.mkdir(self.outputfolder)
            sl(0.5)
        else:
            os.mkdir(self.outputfolder)
            sl(0.5)

    def test_with_dummy_week_investigation(self):

        self.maxDiff = 5000
        cwd = os.getcwd()
        args = argparse.Namespace(variant_data = cwd+"/dodge/tests/inputs/fulltest1_2_allele_ap.txt",
                                    inputtype = "allele",
                                    strainmetadata = cwd+"/dodge/tests/inputs/fulltest1_week_allele_metadata.txt",
                                    outputPrefix = cwd+"/dodge/tests/outputs/fulltest1_out",
                                    distances = cwd+"/dodge/tests/inputs/fulltest1_background_pairwise_distances.txt",
                                    # distances = "apg_testing/input_data/stm/vic_only_5min_2018bg_2021-01-26_2021-02-01_pairwise_distances.txt",
                                    inclusters = cwd+"/dodge/tests/inputs/fulltest1_background_all_clusters.txt",
                                    # inclusters = False,
                                    # inclusters = "/Users/mjohnpayne/Library/CloudStorage/OneDrive-UNSW/Salmonella/2_stage_clustering/manuscript/bioinformatics_submission/revision/Aus2months/v7_static_background_all_clusters.txt",
                                    # inclusters = "/Users/mjohnpayne/Library/CloudStorage/OneDrive-UNSW/Salmonella/2_stage_clustering/manuscript/bioinformatics_submission/revision/Aus2months/v4_background_all_clusters.txt",
                                    # startdate = False,
                                    startdate = "2017-01-15",
                                    enddate = "2017-01-31",
                                    # enddate = "2016-12-31",
                                    background_data = False,
                                    no_cores = 1,
                                    enterobase_data = False,
                                    timesegment = 'week',
                                    timewindow = 28,
                                    dist_limits = "1-5",
                                    max_missmatch = 100,
                                    minsize = 5,
                                    outbreakmethod = "dodge",
                                    static_cutoff = 5,
                                    nonomenclatureinid = False,
                                    isolatecolumn = False,
                                    exclude_time_in_static = False,
                                    usegenomes = False)

        dodge.main(args)

        investclusters = args.outputPrefix + "_2017-01-29_2017-02-04_investigation_clusters.txt"

        expected_investclusters = cwd+"/dodge/tests/expected_outputs/fulltest1_expected_investigation_clusters.txt"


        self.assertListEqual(
            list(open(investclusters)),
            list(open(expected_investclusters)))

    def test_with_dummy_month_investigation(self):
        self.maxDiff = 5000
        cwd = os.getcwd()
        args = argparse.Namespace(variant_data = cwd+"/dodge/tests/inputs/fulltest1_2_allele_ap.txt",
                                    inputtype = "allele",
                                    strainmetadata = cwd+"/dodge/tests/inputs/fulltest2_month_allele_metadata.txt",
                                    outputPrefix = cwd+"/dodge/tests/outputs/fulltest2_out",
                                    distances = cwd+"/dodge/tests/inputs/fulltest2_background_pairwise_distances.txt",
                                    # distances = "apg_testing/input_data/stm/vic_only_5min_2018bg_2021-01-26_2021-02-01_pairwise_distances.txt",
                                    inclusters = cwd+"/dodge/tests/inputs/fulltest2_background_all_clusters.txt",
                                    # inclusters = False,
                                    # inclusters = "/Users/mjohnpayne/Library/CloudStorage/OneDrive-UNSW/Salmonella/2_stage_clustering/manuscript/bioinformatics_submission/revision/Aus2months/v7_static_background_all_clusters.txt",
                                    # inclusters = "/Users/mjohnpayne/Library/CloudStorage/OneDrive-UNSW/Salmonella/2_stage_clustering/manuscript/bioinformatics_submission/revision/Aus2months/v4_background_all_clusters.txt",
                                    # startdate = False,
                                    startdate = "2017-03",
                                    enddate = "2017-05",
                                    # enddate = "2016-12-31",
                                    background_data = False,
                                    no_cores = 1,
                                    enterobase_data = False,
                                    timesegment = 'month',
                                    timewindow = 28,
                                    dist_limits = "1-5",
                                    max_missmatch = 100,
                                    minsize = 5,
                                    outbreakmethod = "dodge",
                                    static_cutoff = 5,
                                    nonomenclatureinid = False,
                                    isolatecolumn = False,
                                    exclude_time_in_static = False,
                                    usegenomes = False)

        dodge.main(args)

        investclusters = args.outputPrefix + "_2017-05_investigation_clusters.txt"

        expected_investclusters = cwd+"/dodge/tests/expected_outputs/fulltest2_expected_investigation_clusters.txt"


        self.assertListEqual(
            list(open(investclusters)),
            list(open(expected_investclusters)))

    def tearDown(self):
        shutil.rmtree(self.outputfolder)


if __name__ == '__main__':
    unittest.main()
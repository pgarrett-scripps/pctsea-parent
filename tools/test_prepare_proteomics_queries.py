import argparse
import csv
import tempfile
import unittest
from pathlib import Path

from prepare_proteomics_queries import prepare


class PrepareProteomicsQueriesTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.fasta = self.root / "proteins.fasta"
        self.input_dir = self.root / "input"
        self.output_dir = self.root / "output"
        self.input_dir.mkdir()
        self.fasta.write_text(
            ">sp|P1|ONE_HUMAN Protein one OS=Homo sapiens GN=GENEA\nAAAA\n"
            ">sp|P2|TWO_HUMAN Protein two OS=Homo sapiens GN=GENEA\nAA\n"
            ">sp|P3|THREE_HUMAN Protein three OS=Homo sapiens GN=GENEB\nAAAA\n",
            encoding="utf-8",
        )
        (self.input_dir / "sample.tsv").write_text(
            "Locus\tSpectrum Count\n"
            "sp|P1|ONE_HUMAN\t8\n"
            "sp|P2|TWO_HUMAN\t2\n"
            "sp|P3|THREE_HUMAN\t4\n"
            "contaminant_KERATIN\t9\n"
            "Reverse_sp|P1|ONE_HUMAN\t5\n"
            "sp|M1|ONE_MOUSE\t7\n"
            "sp|MISSING|MISSING_HUMAN\t3\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def args(self, aggregation: str = "max") -> argparse.Namespace:
        return argparse.Namespace(
            input_dir=self.input_dir,
            fasta=self.fasta,
            output_dir=self.output_dir,
            target_suffix="HUMAN",
            aggregation=aggregation,
            nsaf_scale=1_000_000.0,
            top_genes=None,
        )

    def read_query(self, kind: str) -> dict[str, float]:
        with (self.output_dir / kind / "sample.tsv").open(encoding="utf-8") as handle:
            return {
                row["gene"]: float(row["value"])
                for row in csv.DictReader(handle, delimiter="\t")
            }

    def test_max_aggregation_and_nsaf_length_correction(self) -> None:
        summaries = prepare(self.args())
        self.assertEqual(self.read_query("raw"), {"GENEA": 8.0, "GENEB": 4.0})
        nsaf = self.read_query("nsaf")
        self.assertAlmostEqual(nsaf["GENEA"], 2 / 3 * 1_000_000, places=5)
        self.assertAlmostEqual(nsaf["GENEB"], 1 / 3 * 1_000_000, places=5)
        summary = summaries[0]
        self.assertEqual(summary.input_rows, 7)
        self.assertEqual(summary.mapped_rows, 3)
        self.assertEqual(summary.output_genes, 2)
        self.assertEqual(summary.collapsed_rows, 1)
        self.assertEqual(summary.contaminant_rows, 1)
        self.assertEqual(summary.decoy_rows, 1)
        self.assertEqual(summary.foreign_species_rows, 1)
        self.assertEqual(summary.missing_accession_rows, 1)
        self.assertTrue((self.output_dir / "provenance.tsv").is_file())

    def test_sum_aggregation_is_explicit(self) -> None:
        prepare(self.args(aggregation="sum"))
        self.assertEqual(self.read_query("raw"), {"GENEA": 10.0, "GENEB": 4.0})
        nsaf = self.read_query("nsaf")
        self.assertAlmostEqual(nsaf["GENEA"], 750_000.0)
        self.assertAlmostEqual(nsaf["GENEB"], 250_000.0)

    def test_top_genes_limits_each_output(self) -> None:
        args = self.args()
        args.top_genes = 1
        prepare(args)
        self.assertEqual(self.read_query("raw"), {"GENEA": 8.0})
        self.assertEqual(set(self.read_query("nsaf")), {"GENEA"})


if __name__ == "__main__":
    unittest.main()

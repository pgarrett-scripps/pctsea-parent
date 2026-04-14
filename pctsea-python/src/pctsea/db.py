"""MongoDB database layer for PCTSEA."""

from __future__ import annotations

import logging
from typing import Any

from pymongo import MongoClient
from pymongo.cursor import Cursor

log = logging.getLogger(__name__)


class PctseaDB:
    """Interface to the PCTSEA MongoDB database."""

    def __init__(
        self,
        host: str = "localhost",
        port: int = 27017,
        db_name: str = "single_cells_db",
    ) -> None:
        self.client: MongoClient[dict[str, Any]] = MongoClient(host, port)
        self.db = self.client[db_name]
        self.expressions = self.db["expression"]
        self.single_cells = self.db["singleCell"]
        self.datasets = self.db["dataset"]
        self.cell_type_and_gene = self.db["cellTypeAndGene"]

    def close(self) -> None:
        self.client.close()

    # --- Dataset queries ---

    def find_dataset_by_tag(self, tag: str) -> dict[str, Any] | None:
        return self.datasets.find_one({"tag": tag})

    def find_dataset_by_name(self, name: str) -> dict[str, Any] | None:
        return self.datasets.find_one({"name": name})

    def get_all_datasets(self) -> list[dict[str, Any]]:
        return list(self.datasets.find())

    # --- SingleCell queries ---

    def get_single_cells(self, dataset_tags: list[str] | None = None) -> Cursor[dict[str, Any]]:
        query: dict[str, Any] = {}
        if dataset_tags:
            query["datasetTag"] = {"$in": dataset_tags}
        return self.single_cells.find(query).batch_size(10000)

    def count_cells_by_type(self, cell_type: str) -> int:
        return self.single_cells.count_documents({"type": cell_type})

    def count_cells_by_dataset(self, dataset_tag: str) -> int:
        return self.single_cells.count_documents({"datasetTag": dataset_tag})

    # --- Expression queries ---

    def get_expressions_by_gene(
        self, gene: str, dataset_tags: list[str] | None = None
    ) -> Cursor[dict[str, Any]]:
        query: dict[str, Any] = {"gene": gene}
        if dataset_tags:
            query["projectTag"] = {"$in": dataset_tags}
        return self.expressions.find(query).batch_size(10000)

    def get_expressions_by_genes(
        self, genes: list[str], dataset_tags: list[str] | None = None
    ) -> Cursor[dict[str, Any]]:
        query: dict[str, Any] = {"gene": {"$in": genes}}
        if dataset_tags:
            query["projectTag"] = {"$in": dataset_tags}
        return self.expressions.find(query).batch_size(10000)

    def count_expressions_by_gene(self, gene: str, dataset_tags: list[str] | None = None) -> int:
        query: dict[str, Any] = {"gene": gene}
        if dataset_tags:
            query["projectTag"] = {"$in": dataset_tags}
        return self.expressions.count_documents(query)

    def count_expressions_by_gene_celltype_dataset(
        self, gene: str, cell_type: str, dataset_tag: str
    ) -> int:
        return self.expressions.count_documents(
            {
                "gene": gene,
                "cellType": cell_type,
                "projectTag": dataset_tag,
            }
        )

    # --- CellTypeAndGene queries ---

    def find_cell_type_gene(
        self, dataset_tag: str, cell_type: str, gene: str
    ) -> dict[str, Any] | None:
        return self.cell_type_and_gene.find_one(
            {
                "datasetTag": dataset_tag,
                "cellType": cell_type,
                "gene": gene,
            }
        )

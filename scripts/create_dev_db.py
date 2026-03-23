"""Create an empty equilibrator-compatible SQLite database for development.

Usage:
    python scripts/create_dev_db.py                          # default output path
    python scripts/create_dev_db.py --output path/to/dev.sqlite
    python scripts/create_dev_db.py --source path/to/prod.sqlite  # custom source for registries

The script copies the schema from equilibrator_cache's SQLAlchemy models and
seeds the registries table from the production database so that
LocalCompoundCache can load the new file without issues.
"""

import argparse
import sqlite3
from pathlib import Path

from sqlalchemy import create_engine

from equilibrator_cache.models.compound import Base

DEFAULT_SOURCE = Path("artifacts/.cache/equilibrator/compounds.sqlite")
DEFAULT_OUTPUT = Path("artifacts/.cache/equilibrator/dev_compounds.sqlite")


def create_dev_db(source: Path, output: Path) -> None:
    if output.exists():
        output.unlink()
        print(f"Removed existing {output}")

    # Create empty database with full schema
    dev_engine = create_engine(f"sqlite:///{output}", future=True)
    Base.metadata.create_all(dev_engine)
    dev_engine.dispose()
    print(f"Created schema in {output}")

    # Copy registries via raw SQL to avoid ORM validation issues
    src_conn = sqlite3.connect(source)
    dst_conn = sqlite3.connect(output)
    rows = src_conn.execute("SELECT * FROM registries").fetchall()
    cols = [d[0] for d in src_conn.execute("SELECT * FROM registries LIMIT 0").description]
    placeholders = ", ".join(["?"] * len(cols))
    dst_conn.executemany(f"INSERT INTO registries ({', '.join(cols)}) VALUES ({placeholders})", rows)
    dst_conn.commit()
    src_conn.close()
    dst_conn.close()

    print(f"Seeded {len(rows)} registries from {source}")
    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path, default=DEFAULT_SOURCE,
        help="Production .sqlite to copy registries from",
    )
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT,
        help="Output path for the new dev database",
    )
    args = parser.parse_args()
    create_dev_db(args.source, args.output)

__author__ = "alvaro barbeira"
import logging
import sqlite3
def model_indexes(conn):
    cursor = conn.cursor()
    cursor.execute("CREATE INDEX gene_model_summary ON extra (gene)")
    cursor.execute("CREATE INDEX weights_gene ON weights (gene)")
    cursor.execute("CREATE INDEX weights_rsid ON weights (rsid)")
    cursor.execute("CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")

def create_model_db_(conn, extra, weights):
    logging.log(9, "Saving extra")
    extra.to_sql("extra", conn, index=False)
    logging.log(9, "Saving weights")
    weights.to_sql("weights", conn, index=False)
    logging.log(9, "Creating indexes")
    model_indexes(conn)

def create_model_db(path, extra, weights):
    with sqlite3.connect(path) as conn:
        logging.log(9, "Saving extra")
        extra.to_sql("extra", conn, index=False)
        logging.log(9, "Saving weights")
        weights.to_sql("weights", conn, index=False)
        logging.log(9, "Creating indexes")
        model_indexes(conn)
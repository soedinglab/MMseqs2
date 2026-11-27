import pandas as pd
import argparse


def delta(col1, col2, prec):
    df = pd.DataFrame()
    df['diff'] = abs(col1 - col2)
    return df

def test_differences(parquet_file, tsv_file):

    parquet = pd.read_parquet(parquet_file)
    tsv = pd.read_csv(tsv_file, sep = '\t', header = 0)
    different = False

    # Blob to String
    for name in parquet.select_dtypes(include='object').columns.tolist():
        parquet[name] = parquet[name].str.decode('utf-8')

    parquet = parquet.sort_values(by = ['query','target'],ignore_index = True)
    tsv = tsv.sort_values(by = ['query','target'],ignore_index = True)

    # Check row numbers
    if (parquet.size != tsv.size):
        print("Different table size")
        print(parquet.size(), tsv.size())
        different = True


    # Check for NaNs difference
    if (parquet.isna().values.any() ^ tsv.isna().values.any()):
        print("Dissagrement on existance of NaN values")
        print(parquet.isna().any())
        print(tsv.isna().any())
        different = True

    # Compare values
    for p,t in zip(parquet.select_dtypes(exclude='object').columns.tolist(), tsv.select_dtypes(exclude='object').columns.tolist()):
        prec = 0.001
        dff = delta(parquet[p],tsv[t],prec)

        if sum(dff['diff'] >= prec) > 0:
            parq = parquet[dff['diff'] >= prec][p]
            tsvs = tsv[dff['diff'] >= prec][t]
            pout = pd.concat([parquet[dff['diff'] >= prec]["query"], parquet[dff['diff'] >= prec]["target"], parq,tsvs], axis = 1)
            pout.columns=["query","target",f"Parquet {p}", f"TSV {t}"]
            print(f"Difference bigger than {prec} in column {p}")
            print(pout.head())
            different = True

    # Compare Strings
    for p,t in zip(parquet.select_dtypes(include='object').columns.tolist(), tsv.select_dtypes(include='object').columns.tolist()):
        diff = parquet[p] != tsv[t]
        if (sum(diff) > 0):
            parq = parquet[diff][p]
            tsvs = tsv[diff][t]
            pout = pd.concat([parq,tsvs], axis = 1)
            pout.columns=[f"Parquet {p}", f"TSV {t}"]
            print(pout.head())
            different = True

    print(parquet)
    print(tsv)
    
    if different == False:
        print("No difference in table values found")
        return 1
    else:
        return 0


def main():
    parser = argparse.ArgumentParser(description="Checks for table differences")
    parser.add_argument("parquet", help="Parquet file", type=str)  # Positional
    parser.add_argument("tsv", help="TSV file", type=str)
    args = parser.parse_args()
    return test_differences(args.parquet, args.tsv)

main()

# import argparse
# import duckdb
# PREC = 0.001  # numeric tolerance
# def test_differences(parquet_file, tsv_file):
#     con = duckdb.connect()
#     # Read both tables lazily from disk
#     con.execute(f"""
#         CREATE TEMP TABLE parquet_tbl AS
#         SELECT * FROM read_parquet('{parquet_file}');
#     """)
#     con.sql("SELECT * FROM parquet_tbl").show()
#     con.execute(f"""
#         CREATE TEMP TABLE tsv_tbl AS
#         SELECT * FROM read_csv_auto('{tsv_file}',
#                                     delim='\t',
#                                     header=True);
#     """)
#     con.sql("SELECT * FROM tsv_tbl").show()
#     different = False
#     # Get column info (assume same column names in both)
#     cols = con.execute("""
#         SELECT p.column_name,
#                p.data_type  AS parquet_type,
#                t.data_type  AS tsv_type
#         FROM information_schema.columns p
#         JOIN information_schema.columns t
#           ON p.column_name = t.column_name
#          AND p.table_schema = 'temp'
#          AND t.table_schema = 'temp'
#         WHERE p.table_name = 'parquet_tbl'
#           AND t.table_name = 'tsv_tbl'
#         ORDER BY p.ordinal_position
#     """).fetchall()
#     col_names = [c[0] for c in cols]
#     # Check table "size" (number of rows and columns)
#     parquet_rows, parquet_cols = con.execute(
#         "SELECT COUNT(*), COUNT(*) FROM parquet_tbl, (SELECT COUNT(*) FROM information_schema.columns WHERE table_schema='temp' AND table_name='parquet_tbl') c"
#     ).fetchone()
#     tsv_rows, tsv_cols = con.execute(
#         "SELECT COUNT(*), COUNT(*) FROM tsv_tbl, (SELECT COUNT(*) FROM information_schema.columns WHERE table_schema='temp' AND table_name='tsv_tbl') c"
#     ).fetchone()
#     if parquet_rows != tsv_rows or parquet_cols != tsv_cols:
#         print("Different table size")
#         print("Parquet rows/cols:", parquet_rows, parquet_cols)
#         print("TSV rows/cols    :", tsv_rows, tsv_cols)
#         different = True
#     # Build OR-expression for "any NULL" per table (approx of pandas .isna().values.any())
#     any_null_expr_parquet = " OR ".join([f'"{c}" IS NULL' for c in col_names]) or "FALSE"
#     any_null_expr_tsv = " OR ".join([f'"{c}" IS NULL' for c in col_names]) or "FALSE"
#     pq_has_null, tsv_has_null = con.execute(f"""
#         SELECT
#           EXISTS (SELECT 1 FROM parquet_tbl WHERE {any_null_expr_parquet}) AS pq_has_null,
#           EXISTS (SELECT 1 FROM tsv_tbl     WHERE {any_null_expr_tsv})     AS tsv_has_null;
#     """).fetchone()
#     if bool(pq_has_null) ^ bool(tsv_has_null):
#         print("Disagreement on existence of NaN/NULL values")
#         print("parquet has NULLs:", pq_has_null)
#         print("tsv has NULLs    :", tsv_has_null)
#         different = True
#     # Separate numeric vs string-like columns, ignoring blobs here (your TSV has no blobs)
#     numeric_cols = []
#     string_cols = []
#     for name, ptype, ttype in cols:
#         pt = ptype.upper()
#         tt = ttype.upper()
#         if any(x in pt for x in ["DOUBLE", "FLOAT", "DECIMAL", "NUMERIC", "REAL"]) or \
#            any(x in tt for x in ["DOUBLE", "FLOAT", "DECIMAL", "NUMERIC", "REAL"]):
#             numeric_cols.append(name)
#         elif "CHAR" in pt or "VARCHAR" in pt or "STRING" in pt or \
#              "CHAR" in tt or "VARCHAR" in tt or "STRING" in tt:
#             string_cols.append(name)
#         else:
#             # Treat other non-int numerics as numeric; ints will work fine in numeric comparison
#             if "INT" in pt or "INT" in tt:
#                 numeric_cols.append(name)
#             else:
#                 string_cols.append(name)
#     # We assume 'query' and 'target' exist and are join keys as in your original code
#     key_cols = ['query', 'target']
#     # Numeric comparison (difference >= PREC)
#     for col in numeric_cols:
#         if col in key_cols:
#             continue
#         # join on query,target to align rows
#         res = con.execute(f"""
#             SELECT
#               COUNT(*) AS cnt
#             FROM (
#               SELECT
#                 ABS(p."{col}" - t."{col}") AS diff
#               FROM parquet_tbl p
#               JOIN tsv_tbl t
#                 ON p.query = t.query
#                AND p.target = t.target
#             ) d
#             WHERE diff >= {PREC}
#         """).fetchone()
#         cnt = res[0]
#         if cnt and cnt > 0:
#             print(f"Difference bigger than {PREC} in column {col}")
#             # show some example rows
#             sample = con.execute(f"""
#                 SELECT
#                   p.query,
#                   p.target,
#                   p."{col}" AS parquet_value,
#                   t."{col}" AS tsv_value
#                 FROM parquet_tbl p
#                 JOIN tsv_tbl t
#                   ON p.query = t.query
#                  AND p.target = t.target
#                 WHERE ABS(p."{col}" - t."{col}") >= {PREC}
#                 LIMIT 5;
#             """).fetchall()
#             for row in sample:
#                 print(row)
#             different = True
#     # String comparison
#     for col in string_cols:
#         if col in key_cols:
#             continue
#         res = con.execute(f"""
#             SELECT COUNT(*) AS cnt
#             FROM parquet_tbl p
#             JOIN tsv_tbl t
#               ON p.query = t.query
#              AND p.target = t.target
#             WHERE p."{col}" IS DISTINCT FROM t."{col}";
#         """).fetchone()
#         cnt = res[0]
#         if cnt and cnt > 0:
#             print(f"String differences in column {col}")
#             sample = con.execute(f"""
#                 SELECT
#                   p.query,
#                   p.target,
#                   p."{col}" AS parquet_value,
#                   t."{col}" AS tsv_value
#                 FROM parquet_tbl p
#                 JOIN tsv_tbl t
#                   ON p.query = t.query
#                  AND p.target = t.target
#                 WHERE p."{col}" IS DISTINCT FROM t."{col}"
#                 LIMIT 5;
#             """).fetchall()
#             for row in sample:
#                 print(row)
#             different = True
#     if not different:
#         print("No difference in table values found")
#         return 1
#     else:
#         return 0
# def main():
#     parser = argparse.ArgumentParser(description="Checks for table differences using DuckDB")
#     parser.add_argument("parquet", help="Parquet file", type=str)
#     parser.add_argument("tsv", help="TSV file", type=str)
#     args = parser.parse_args()
#     return test_differences(args.parquet, args.tsv)
# if __name__ == "__main__":
#     main()
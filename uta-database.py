# import psycopg2
# import psycopg2.extras
# conn = psycopg2.connect(
#     "host=uta.biocommons.org dbname=uta user=anonymous password=anonymous")
# cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
# cur.execute("select * from uta_20140210.tx_def_summary_v where hgnc='BRCA1'")
# row = cur.fetchone()
# dict(row)
# {
#     'hgnc': 'BRCA1',
#     'tx_ac': 'ENST00000586385',
#     'cds_md5': '5d405c70b9b79add38d28e5011a6ddc0',
#     'es_fingerprint': '95d60b8d62f5c23cbeff3499eedf5e89',
#     'cds_start_i': 144,
#     'cds_end_i': 666,
#     'starts_i': [0, 148, 226, 267, 351, 406, 480, 541],
#     'ends_i': [148, 226, 267, 351, 406, 480, 541, 781],
#     'lengths': [148, 78, 41, 84, 55, 74, 61, 240],
# }

import psycopg2

conn = psycopg2.connect(
    "host=uta.biocommons.org dbname=uta user=anonymous password=anonymous")

cur = conn.cursor()

# Get the names of all tables in the public schema
cur.execute("""
    SELECT table_name
    FROM information_schema.tables
""")
rows = cur.fetchall()

# Print the names of all tables
for row in rows:
    print(row[0])
# Append all these names into a txt file line by line
with open('tables.txt', 'w') as f:
    for row in rows:
        f.write(row[0] + '\n')

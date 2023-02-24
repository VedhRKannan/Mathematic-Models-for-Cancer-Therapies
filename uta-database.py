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

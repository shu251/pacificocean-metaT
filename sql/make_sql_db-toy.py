import sqlite3
import pandas as pd

conn = sqlite3.connect('/scratch/user/skhu/SPOT-ALOHA/toy-metat.db')

tpm_data = pd.read_csv('/scratch/user/skhu/SPOT-ALOHA/03-abundance_tables/TPMTable_ByContig_150k.csv')


## Set tpm_data object to SQL:
# tpm_data.to_sql, table name = "TPM", conn = sqlite db, if exists replace and do not index

tpm_data.to_sql('TPM',conn, if_exists='replace',index=False)


# Add metadata

metadata = pd.read_csv('/home/skhu/pacificocean-metaT/input-data/complete-sample-list.txt')
metadata.to_sql('metadata', conn, if_exists='replace', index=False)


# Add annotation data
taxfxn = pd.read_csv('/scratch/user/skhu/SPOT-ALOHA/02-annotation_table/TaxonomicAndFunctionalAnnotations.csv')
taxfxn.to_sql('annot', conn, if_exists='replace', index=False)

conn.close()


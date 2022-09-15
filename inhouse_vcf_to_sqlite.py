import argparse
import pandas as pd
import sqlite3
import gzip
import sys
import csv
import re
import os


class Vcf_to_sqlite:
    def __init__(self):
        self.sample_vcf = args.vcf
        self.db = 'inhouse.db'
        self.sample_name = args.sample
        self.old_inhouse = args.import_inhouse_vcf
        self.iteratefolder = args.add_folder

    def create_empty_db(self):
        # CREATE DB IF NOT EXISTS
        conn = sqlite3.connect('inhouse.db')
        print "Opened database successfully"

        conn.execute('''CREATE TABLE IF NOT EXISTS VARIANTS
                (ID TEXT PRIMARY KEY,
                CHROM     CHAR(5)    NOT NULL,
                POS       INT     NOT NULL,
                var_ID    CHAR(50),
                REF       CHAR(50),
                ALT       CHAR(50),
                QUAL      REAL);''')
        print "Table created successfully"

        conn.execute('''CREATE TABLE IF NOT EXISTS SAMPLES
                (ID INTEGER PRIMARY KEY AUTOINCREMENT,
                VARIANT TEXT,
                FILTER  TEXT,
                SAMPLE  TEXT,
                GT      CHAR(50),
                UNIQUE(VARIANT, SAMPLE));''')

        print "Table created successfully"

        conn.close()
        print "Database closed successfully"

    def import_sample_vcf(self, infile, sample):
        types_dict = {0: str, 1: int, 2:str , 3: str, 4: str, 5: str, 6: str}
        vcf = pd.read_csv(open(infile,'rb'), index_col=False,delimiter='\t',comment='#', header=0, usecols=[0, 1, 2, 3, 4, 5, 6, 9], dtype=types_dict ,names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"], converters = { 'INFO' : lambda x : x.split(':')[0] })
        conn = sqlite3.connect(self.db)
        print sample
        print "Opened database successfully"
        cursor = conn.cursor()
        for index, row in vcf.iterrows():
            # if not exist
            insert_with_param = """INSERT OR IGNORE INTO VARIANTS (ID, CHROM, POS, var_ID, REF, ALT, QUAL) VALUES (?,?,?,?,?,?,?)"""
            data_tuple = ("{}{}{}{}".format(row['CHROM'],row['POS'],row['REF'],row['ALT']), row['CHROM'], row['POS'], row['ID'], row['REF'], row['ALT'], row['QUAL'])
            cursor.execute(insert_with_param, data_tuple)

            # Update sample data:
            insert_sample_with_param = """INSERT INTO SAMPLES (VARIANT, FILTER, GT, SAMPLE) VALUES (?, ?, ?, ?)"""
            sample_tuple = ("{}{}{}{}".format(row['CHROM'],row['POS'],row['REF'],row['ALT']), row['FILTER'], row['INFO'], sample)
            cursor.execute(insert_sample_with_param, sample_tuple)
            conn.commit()

        print "Records created successfully"
        cursor.close()
        conn.close()


    def printstats(self):
        conn = sqlite3.connect(self.db)
        print "Opened database successfully:"
        cursor = conn.cursor()
        # Count all variants
        cursor.execute('SELECT COUNT(*) from VARIANTS')
        cur_result = cursor.fetchone()
        print("There are a total of {} variants in the db.".format(cur_result[0]))
        # Count all samples
        cursor.execute('SELECT COUNT(DISTINCT(SAMPLE)) from SAMPLES')
        cur_result = cursor.fetchone()
        print("There are a total of {} samples in the db.".format(cur_result[0]))



    def print_vcf(self):
        conn = sqlite3.connect(self.db)
        query = """SELECT ID,CHROM,POS,var_ID,REF,ALT,QUAL,SAMPS,HETS,HOMS FROM VARIANTS
        LEFT JOIN ( 
            SELECT VARIANT,GROUP_CONCAT(SAMPLE||'('||CASE WHEN GT == '0/1' THEN 'HET' ELSE 'HOM' END||')','_') AS SAMPS,GT,SUM(CASE WHEN GT == '0/1' THEN 1 ELSE 0 END) AS HETS, SUM(CASE WHEN GT == '1/1' THEN 1 ELSE 0 END) AS HOMS FROM SAMPLES GROUP BY VARIANT
        ) ON VARIANTS.ID = VARIANT"""
        df = pd.read_sql_query(query, conn)
        writer = csv.writer(open('hg19_inhouse.vcf','wb'), delimiter="\t", quoting=csv.QUOTE_NONE, quotechar="", escapechar=' ')
        df[['HETS', 'HOMS']] = df[['HETS','HOMS']].fillna(value=0)
        df[['QUAL']] = df[['QUAL']].fillna(value=".")
        header = ['##fileformat=VCFv4.1',
        '##fileDate=20170821',
        '##source=CentoMD',
        '##reference=GRChr37',
        '##INFO=<ID=IHSAMPLES,Number=1,Type=String,Description="Sample names from inhouse db that contains this variant and het state in parentheses">',
        '##INFO=<ID=IHHET,Number=1,Type=String,Description="Number of inhouse hets for this variant">',
        '##INFO=<ID=IHHOM,Number=1,Type=String,Description="Number of inhouse homozygotes for this variant">']
        
        lasthead = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
        
        
        
        for i in header:
            writer.writerow([i])
        writer.writerow(lasthead)
        for index, row in df.iterrows():
            if row[5] != "0":
                writestring = [row[1],row[2], row[3], row[4], row[5], row[6], "PASS", "IHSAMPLES={};IHHET={};IHHOM={}".format( row[7], int(row[8]), int(row[9]))]
                writer.writerow(writestring)

    def get_from_inhouse(self):
        ''' Get variants from previous inhouse vcf into db:'''
        
        conn = sqlite3.connect(self.db)
        print "Opened database successfully"
        cursor = conn.cursor()
        
        # Read gz line by line
        with gzip.open(self.old_inhouse, mode="rt") as f:
            lines = filter(None, (line.rstrip() for line in f))
            for line in lines:
                if not line.startswith("#"):
            
                    # INSERT VARIANTS
                    insert = """INSERT OR IGNORE INTO VARIANTS (ID, CHROM, POS, REF, ALT, var_ID) VALUES (?,?,?,?,?,?)"""
                    data_tuple = ("{}{}{}{}".format(line.split("\t")[0],line.split("\t")[1],line.split("\t")[3],line.split("\t")[4]),line.split("\t")[0],line.split("\t")[1],line.split("\t")[3],line.split("\t")[4], line.split("\t")[2])
                    cursor.execute(insert, data_tuple)
                    conn.commit()
            
                    # INSERT SAMPLES:
                    
                    length = int(line.split("\t")[-1].split(";")[1].split("=")[1]) + int(line.split("\t")[-1].split(";")[2].split("=")[1])
                    
                    for l in range(length):
                        gt = "0/1" if line.split("\t")[-1].split(";")[0].replace("IHSAMPLES=","").split(")_")[l].split("(")[1].replace(")","") == "HET" else "1/1"
                        insert_s = """INSERT OR IGNORE INTO SAMPLES (VARIANT, FILTER, SAMPLE, GT) VALUES (?,?,?,?)"""
                        data_tuple_s = ("{}{}{}{}".format(line.split("\t")[0],line.split("\t")[1],line.split("\t")[3],line.split("\t")[4]),"PASS",
                                    line.split("\t")[-1].split(";")[0].replace("IHSAMPLES=","").split(")_")[l].split("(")[0],
                                    gt)
                        
                        cursor.execute(insert_s, data_tuple_s)
                        conn.commit()
        print "Records created successfully"
        cursor.close()
        conn.close()

    def iterate_folder(self):
        samples = {}
        sample = ""
        path = self.iteratefolder
        for path, subdirs, files in os.walk(path):
            for name in files:
                if re.search(".raw.variants.Filtered.individual.annotated.subset.", name):
                    sample = name.split('.')[0]
                    # Run insert into db:
                    samples[sample] = os.path.join(path,name)
        return samples




if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''Takes a vcf, enters into an existing sqlite-db, or if non-existing - creates one 
    You have to sort the output: 
    docker run -it --rm \  
    -v $PWD:/data \  
    -v /illumina/runs_diag/prod_pipeline/genomes/H_sapiens/b37/:/ref \  
    broadinstitute/gatk:4.1.4.1 \  
    gatk SortVcf \  
    --INPUT /data/hg19_inhouse.vcf \  
    --OUTPUT /data/out.vcf \  
    --SEQUENCE_DICTIONARY /ref/human_g1k_v37.dict
    
    Remember to bgzip & tabix 
    
    ''', epilog='')
    parser.add_argument('-v', '--vcf',                  help="Name of pr sample vcf to be imported into database",                                                      required=False                          )
    parser.add_argument('-c', '--create_db',            help='Creates an empty sqlite-db',                                                                              required=False,     action='store_true' )
    parser.add_argument('-s', '--sample',               help="Sample name of the vcf",                                                                                  required=False                          )
    parser.add_argument('-t', '--stats',                help="Print stats",                                                                                             required=False,     action='store_true' )    
    parser.add_argument('-p', '--print_vcf',            help="output vcf for use in pipeline",                                                                          required=False,     action='store_true' )
    parser.add_argument('-i', '--import_inhouse_vcf',   help="Import a vcf with inhouse variants",                                                                      required=False                          )
    parser.add_argument('-f', '--add_folder',           help="Add a folder, and search it for files containing *.raw.variants.Filtered.individual.annotated.subset.* ", required=False                          )
    

    args = parser.parse_args()
    run = Vcf_to_sqlite()
    if args.create_db:
        run.create_empty_db()
        sys.exit()
    if args.stats:
        run.printstats()
        sys.exit()
    if args.print_vcf:
        run.print_vcf()
        sys.exit()
    if args.import_inhouse_vcf:
        run.get_from_inhouse()
        sys.exit()
    if args.add_folder:
        samples = run.iterate_folder()
        for sample_name, sample_path in samples.items():
            print sample_name
            print sample_path
            run.import_sample_vcf(sample_path, sample_name)
        sys.exit()
    if args.vcf:
        run.import_sample_vcf(run.sample_vcf, run.sample_name)

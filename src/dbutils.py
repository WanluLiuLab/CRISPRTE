# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#
# ------------------------- #
# Python Modules
# ------------------------- #

import sqlite3
import json
from abc import ABC
import psycopg2
# rom .utils import SymbolToNumber, NumberToSymbol, PatternToNumber, NumberToPattern

class SQLBase(ABC):

    @staticmethod
    def connect(db, **args):
        return sqlite3.connect(db, **args)

    @staticmethod
    def delete_gtf_table(conn,tbname):
        c = conn.cursor()
        c.execute("""DROP TABLE {} """.format(tbname))
        conn.commit()
        conn.close()

    @staticmethod
    def create_ttf_table(conn):
        ttf_format = """ (name BYTEA NOT NULL,
                          copies INTEGER NOT NULL,
                          source BYTEA NOT NULL,
                          tid INTEGER NOT NULL,
                          dfamId INTEGER NOT NULL,
                          classification BYTEA NOT NULL,
                          clades BYTEA NOT NULL,
                          description BYTEA NOT NULL,
                          length INTEGER NOT NULL,
                          consensus BYTEA NOT NULL)
                     """
        cursor = conn.cursor()
        cursor.execute("""CREATE TABLE ttf """ + ttf_format)
        conn.commit()
        conn.close()

    @staticmethod
    def create_dtf_table(conn):
        dtf_format = """ (tid INTEGER NOT NULL,
                          dup INTEGER NOT NULL,
                          optimal_alignment_score INTEGER NOT NULL,
                          query_begin INTEGER NOT NULL,
                          query_end INTEGER NOT NULL,
                          target_begin INTEGER NOT NULL,
                          target_end INTEGER NOT NULL,
                          target_sequence BYTEA NOT NULL,
                          cigar BYTEA NOT NULL)
                     """
        cursor = conn.cursor()
        cursor.execute("""CREATE TABLE dtf """ + dtf_format)
        conn.commit()
        conn.close()

    @staticmethod
    def create_dtf_index(conn):
        sqls = ["CREATE INDEX tid on dtf(tid)", "CREATE INDEX dup on dtf(dup)"]
        c = conn.cursor()
        for sql in sqls:
            c.execute(sql)
        conn.commit()
        conn.close()

    @staticmethod
    def create_gtf_table(conn, tbname):
        gtf_format = """ (seqname BYTEA NOT NULL,
                        source BYTEA NOT NULL,
                        feature BYTEA NOT NULL,
                        starting INTEGER NOT NULL,
                        ending INTEGER NOT NULL,
                        score BYTEA NOT NULL,
                        strand VARCHAR(1) NOT NULL,
                        frame VARCHAR(1) NOT NULL,
                        gene_id BYTEA,
                        gene_version BYTEA,
                        transcript_id BYTEA,
                        transcript_version BYTEA,
                        transcript_name BYTEA,
                        transcript_source BYTEA,
                        transcript_biotype BYTEA,
                        transcript_support_level BYTEA,
                        exon_number BYTEA,
                        exon_id BYTEA,
                        exon_version BYTEA,
                        exon_name BYTEA,
                        exon_source BYTEA,
                        exon_biotype BYTEA,
                        protein_id BYTEA,
                        protein_version BYTEA,
                        protein_name BYTEA,
                        protein_source BYTEA,
                        protein_biotype BYTEA,
                        ccds_id BYTEA,
                        ccds_version BYTEA,
                        ccds_name BYTEA,
                        ccds_source BYTEA,
                        ccds_biotype BYTEA,
                        gene_name BYTEA,
                        gene_source BYTEA,
                        gene_biotype BYTEA,
                        family_id BYTEA,
                        class_id BYTEA,
                        tag BYTEA)"""
        cursor = conn.cursor()
        cursor.execute("""CREATE TABLE {} """.format(tbname) + gtf_format)
        conn.commit()
        conn.close()


    @staticmethod
    def delete_table(conn, tbname):
        c = conn.cursor()
        c.execute("""DROP TABLE {} """.format(tbname))
        conn.commit()
        conn.close()

    @staticmethod
    def create_table1(conn):
        pass

    @staticmethod
    def create_table1_index(conn):
        sqls = ["CREATE INDEX GidIndex on table1(gid)", "CREATE INDEX ClassIndex on table1(te_class)","CREATE INDEX DupIndex on table1(te_dup)"]
        c = conn.cursor()
        for sql in sqls:
            c.execute(sql)
        conn.commit()
        conn.close()

    @staticmethod
    def create_table2(conn):
        pass
        
    @staticmethod
    def insert_gtf_table(cursor, tabname, data):
        for i in data:
            const, attr = i[0], i[1]
            sql = "INSERT INTO {}(seqname, source, feature, start, end, score, strand, frame, ".format(tabname) + ','.join(attr.keys()) + ") VALUES('%s', '%s', '%s','%d', %d, '%s', '%s', '%s', " + ",".join(["%s"] * len(attr)) + ")"
            cursor.execute(sql % tuple(const[:3] + [int(const[3]) ,int(const[4])] + const[5:]  + list(map(lambda x:"'{}'".format(x), attr.values()))))

    @staticmethod
    def insert_table1(cursor, data):
        sql = "INSERT INTO table1 (pos, gid, pam, upstream, downstream, gscore_moreno, gscore_azimuth, anno_class, te_class, te_dup) VALUES('%s', %d, '%s','%s', '%s', %f, %f, '%s', '%s', '%s')"
        cursor.execute(sql % data)


    @staticmethod
    def insert_table2(cursor, data):
        sql = "INSERT INTO table2 (gid, gseq,mm1,mm2,mm3) VALUES(%d, '%s', '%s','%s', '%s')"
        cursor.execute(sql % data)

    @staticmethod
    def insert_table2_nomm(cursor, data):
        sql = "INSERT INTO table2 (gid, gseq) VALUES(%d, '%s')"
        cursor.execute(sql % data)

    @staticmethod
    def insert_table2_gid(cursor, gid, data):
        sql = "UPDATE table2 SET mm1='%s',mm2='%s',mm3='%s' WHERE (gid = {})".format(gid)
        cursor.execute(sql % data)

    @staticmethod
    def select_gid_table1_tedup(cursor, te_dup):
        cursor.execute("SELECT gid FROM table1 WHERE te_dup = '{}'".format(te_dup))
        result = list(map(lambda x:x[0],cursor.fetchall()))
        return result

    @staticmethod
    def select_gid_table1_teclass(cursor, te_class):
        cursor.execute("SELECT gid FROM table1 WHERE te_class = '{}'".format(te_class))
        result = list(map(lambda x:x[0],cursor.fetchall()))
        return result

    @staticmethod
    def select_gid_seqinfo_table1_tedup(cursor, te_dup):
        cursor.execute("SELECT gid,pos,gscore_moreno,gscore_azimuth,pam,upstream,downstream FROM table1 WHERE te_dup = '{}'".format(te_dup))
        result = list(cursor.fetchall())
        return result

    @staticmethod
    def select_gid_seqinfo_table1_teclass(cursor, te_class):
        cursor.execute("SELECT gid,pos,gscore_moreno,gscore_azimuth,pam,upstream,downstream,te_dup FROM table1 WHERE te_class = '{}'".format(te_class))
        result = list(cursor.fetchall())
        return result

    @staticmethod
    def select_gid_anno_table1_gid(cursor, gid):
        cursor.execute("SELECT gid, anno_class, te_class, te_dup FROM table1 WHERE gid = {}".format(gid))
        result = list(cursor.fetchall())
        return result

    @staticmethod
    def select_gid_anno_table1_gids(cursor, gids):
        cursor.execute("SELECT gid, anno_class, te_class, te_dup FROM table1 WHERE gid IN ({})".format(','.join(map(str, gids))))
        result = list(cursor.fetchall())
        return result


    @staticmethod
    def select_anno_table1_gid(cursor, gid):
        cursor.execute("SELECT anno_class, te_class, te_dup FROM table1 WHERE gid = {}".format(gid))
        result = list(cursor.fetchall())
        return result

    @staticmethod
    def select_allinfo_table1_gid(cursor, gid):
        cursor.execute("SELECT pos, gid, gscore_moreno,gscore_azimuth, pam, upstream, downstream, anno_class, te_class, te_dup FROM table1 WHERE gid = {}".format(gid))
        result = list(cursor.fetchall())
        return result

    @staticmethod
    def select_allinfo_table1_gids(cursor, gids):
        cursor.execute("SELECT pos, gid, gscore_moreno,gscore_azimuth, pam, upstream, downstream, anno_class, te_class, te_dup FROM table1 WHERE gid IN ({})".format(','.join(gids)))
        result = list(cursor.fetchall())
        return result

    @staticmethod
    def select_table2_gid(cursor, gid):
        cursor.execute(
            "SELECT gid, gseq,mm1,mm2,mm3 FROM table2 WHERE gid = {}".format(gid))
        result = list(cursor.fetchall()[0])
        for i in range(2, 5):
            result[i] = json.loads(result[i]) if result[i] != 'None' and result[i] != None and result[i] != "\n" else None
        return result

    @staticmethod
    def select_gseq_table2_gid(cursor, gid):
        cursor.execute(
            "SELECT gseq FROM table2 WHERE gid = {}".format(gid))
        return cursor.fetchall()[0][0]

    @staticmethod
    def select_gseq_table2_gids(cursor, gids):
        cursor.execute(
            "SELECT gseq FROM table2 WHERE gid IN ({})".format(','.join(gids)))
        return list(map(lambda x:x[0], cursor.fetchall()))


    @staticmethod
    def select_gseq_table2_where_mm(cursor):
        cursor.execute(
            "SELECT gid,mm1,mm2,mm3 FROM table2 WHERE mm1 is not null".format(gid))
        result = cursor.fetchall()
        return result

    @staticmethod
    def select_gtf_by_genename(cursor, tabname, genename):
        cursor.execute("SELECT seqname, source, feature, starting, ending, score, strand, frame, gene_name FROM {} WHERE gene_name='{}'".format(tabname, genename))
        result = cursor.fetchall()
        return list(map(lambda x:list(x)[:8] + [{"gene_name":x[8]}], result))   

    @staticmethod
    def select_gtf_by_ensembl_geneid(cursor, geneid, tabname = "pcg"):
        cursor.execute("SELECT seqname, source, feature, starting, ending, score, strand, frame, gene_name FROM {} WHERE gene_id='{}'".format(tabname, geneid))
        result = cursor.fetchall()
        return list(map(lambda x:list(x)[:8] + [{"gene_name":x[8]}], result))    

    @staticmethod
    def select_gtf_by_tename(cursor, tename, tabname="te"):
        if "dup" in tename:
            cursor.execute("SELECT seqname, source, feature, starting, ending, score, strand, frame, gene_id, transcript_id FROM {} WHERE transcript_id='{}'".format(tabname, tename))
        else:
            cursor.execute("SELECT seqname, source, feature, starting, ending, score, strand, frame, gene_id, transcript_id FROM {} WHERE gene_id='{}'".format(tabname, tename))
        result = cursor.fetchall()
        return list(map(lambda x:list(x)[:8] + [{"gene_name":x[8]}], result))   
    
    @staticmethod
    def select_gtf_by_region(cursor, tabname, chrom, start, end):
        
        cursor.execute("SELECT seqname, source, feature, starting, ending, score, strand, frame, gene_name FROM {} WHERE seqname='{}' AND start>{} AND end<{}".format(tabname, chrom, start, end))
        result = cursor.fetchall()
        return list(map(lambda x:dict(zip(["seqname","source","feature","start","end","score","strand","frame","name"], x)), result))

    @staticmethod
    def select_te_gtf_by_region(cursor, tabname, chrom, start, end):
        
        cursor.execute("SELECT seqname, source, feature, starting, ending, score, strand, frame, gene_id, transcript_id FROM {} WHERE seqname='{}' AND start>{} AND end<{}".format(tabname, chrom, start, end))
        result = cursor.fetchall()
        return list(map(lambda x:dict(zip(["seqname","source","feature","start","end","score","strand","frame","gene_id","name"], x)), result))

class SQLite(SQLBase):
    @staticmethod
    def create_table1(conn):
        format = """ (pos VARCHAR(50) NOT NULL,
                    gid INTEGER NOT NULL,
                    gscore_moreno DOUBLE NOT NULL,
                    gscore_azimuth DOUBLE NOT NULL,
                    pam VARCHAR(3) NOT NULL,
                    upstream VARCHAR(6) NOT NULL,
                    downstream VARCHAR(6) NOT NULL,
                    anno_class VARCHAR(50) NOT NULL,
                    te_class VARCHAR(40),
                    te_dup VARCHAR(40))"""
        c = conn.cursor()
        c.execute("""CREATE TABLE {} """.format("table1") + format)
        conn.commit()
        conn.close()

    @staticmethod
    def create_table2(conn):
        format = """ (gid INTEGER NOT NULL PRIMARY KEY,
    gseq VARCHAR(20) NOT NULL,
    mm1 BLOB,
    mm2 BLOB,
    mm3 BLOB)"""
        c = conn.cursor()
        c.execute("""CREATE TABLE {} """.format("table2") + format)
        conn.commit()
        conn.close()

class PGSQL(SQLBase):
    @staticmethod
    def connect(dbname):
        print("password:", end='')
        password=input()
        return psycopg2.connect("dbname={} user=postgres password={}".format(dbname, password))

    @staticmethod
    def create_table1(conn):
        format = """ (pos VARCHAR(50) NOT NULL,
                    gid INTEGER NOT NULL,
                    gscore_moreno FLOAT NOT NULL,
                    gscore_azimuth FLOAT NOT NULL,
                    pam VARCHAR(3) NOT NULL,
                    upstream VARCHAR(6) NOT NULL,
                    downstream VARCHAR(6) NOT NULL,
                    anno_class BYTEA NOT NULL,
                    te_class INTEGER, 
                    te_dup INTEGER)"""
        c = conn.cursor()
        c.execute("""CREATE TABLE {} """.format("table1") + format)
        conn.commit()
        conn.close()

    @staticmethod
    def create_table2(conn):
        format = """ (gid INTEGER NOT NULL PRIMARY KEY,
                        gseq VARCHAR(20) NOT NULL,
                        mm1 INTEGER[],
                        mm2 INTEGER[],
                        mm3 INTEGER[])"""
        c = conn.cursor()
        c.execute("""CREATE TABLE {} """.format("table2") + format)
        conn.commit()
        conn.close()

    @staticmethod
    def insert_table1(cursor, data):
        sql = "INSERT INTO table1 (pos, gid, gscore_moreno, gscore_azimuth, pam, upstream, downstream,  anno_class, te_class, te_dup) VALUES('%s', %d,  %f, %f, '%s','%s', '%s', '%s', %d, %d)"
        cursor.execute(sql % data)

    @staticmethod
    def insert_table2(cursor, data):
        sql = "INSERT INTO table2 (gid, gseq, mm1, mm2, mm3) VALUES(%d, '%s', '%s', '%s', '%s')"
        cursor.execute(sql % data)

    @staticmethod
    def insert_gtf_table(cursor, tabname, data):
        for i in data:
            const, attr = i[0], i[1]
            sql = "INSERT INTO {}(seqname, source, feature, starting, ending, score, strand, frame, ".format(tabname) + ','.join(attr.keys()) + ") VALUES('%s', '%s', '%s', %d, %d, '%s', '%s', '%s', " + ",".join(["%s"] * len(attr)) + ")"
            cursor.execute(sql % tuple(const[:3] + [int(const[3]) ,int(const[4])] + const[5:]  + list(map(lambda x:"'{}'".format(x), attr.values()))))

    @staticmethod
    def insert_ttf_table(cursor, data):
        sql = "INSERT INTO ttf (name, copies, source, tid, dfamId, classification, clades, description, length, consensus) VALUES ('%s', %d, '%s', %d, %d, '%s', '%s', '%s', %d, '%s')"
        for i in data:
            cursor.execute(sql % i)

    @staticmethod
    def insert_dtf_table(cursor, data):
        sql =  "INSERT INTO dtf (tid, dup, optimal_alignment_score, query_begin, query_end, target_begin, target_end, target_sequence, cigar) VALUES (%d, %d, %d, %d, %d, %d, %d, '%s','%s')"
        for i in data:
            cursor.execute(sql % i)

    @staticmethod
    def select_te_ttf_by_te(cursor, tid):
        cursor.execute("SELECT name, copies, source, dfamId, classification, clades, description, length, consensus FROM ttf WHERE tid = {}".format(tid))
        result = cursor.fetchall()
        result = dict(zip(["name", "copies", "source", "dfamid", "classification", "clades", "description", "length", "consensus"], result))
        for k,v in result.items():
            if type(v) == memoryview:
                result[k] = v.tobytes.decode()
        return result

    @staticmethod
    def select_te_dtf_by_te(cursor, tid, copy = None):
        if copy:
            cursor.execute("SELECT optimal_alignment_score, query_begin, query_end, target_begin, target_end, target_sequence, cigar FROM dtf WHERE tid = {} and dup = {}".format(tid, copy))
            result = cursor.fetchall()
            result = dict(zip(["optimal_alignment_score", "query_begin", "query_end", "target_begin", "target_end", "target_sequence", "cigar"], result))
            for k,v in result.items():
                if type(v) == memoryview:
                    result[k] = v.tobytes().decode()
        else:
            cursor.execute("SELECT optimal_alignment_score, query_begin, query_end, target_begin, target_end, target_sequence, cigar FROM dtf WHERE tid = {}".format(tid))
            result = cursor.fetchall()
            result = list(map(lambda x:dict(zip(["optimal_alignment_score", "query_begin", "query_end", "target_begin", "target_end", "target_sequence", "cigar"], x)), result))
            for i in range(len(result)):
                for k,v in result[i].items():
                    if type(v) == memoryview:
                        result[i][k] = v.tobytes().decode()
        return result
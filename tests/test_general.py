# -*- coding:utf-8 -*-
from snapgene_reader import snapgene_file_to_seqrecord, snapgene_file_to_gbk, snapgene_file_to_dict
from Bio import SeqIO
import os
import sys
reload(sys)
sys.setdefaultencoding('utf8')

TEST_DIR = os.path.join('test_samples')


def test_parse(tmpdir):
    all_files = [f for f in os.listdir(TEST_DIR) if f.endswith('.dna')]
    assert len(all_files)
    for fname in all_files:
        fpath = os.path.join(TEST_DIR, fname)
        record = snapgene_file_to_seqrecord(fpath)
        assert len(record.seq) > 10
        with open(os.path.join(str(tmpdir), fname + '.gb'), 'w') as f:
            SeqIO.write([record, ], f, 'genbank')


def test_convert_gbk(tmpdir):
    all_files = [f for f in os.listdir(TEST_DIR) if f.endswith('.dna')]
    assert len(all_files)
    for fname in all_files:
        infpath = os.path.join(TEST_DIR, fname)
        read_file_object = open(infpath)
        outfpath = os.path.join(tmpdir, fname)
        write_file_object = open(outfpath + ".gbk", "w")
        snapgene_file_to_gbk(read_file_object, write_file_object)
        read_file_object = open(infpath)
        data = snapgene_file_to_dict(fileobject=read_file_object)
        for key in data:
            print key, data[key]


if __name__ == '__main__':
    test_convert_gbk("./")

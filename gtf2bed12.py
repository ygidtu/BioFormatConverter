#!/usr/bin/env python3
# -*- coding:utf-8 -*-

u"""
convert gtf files to bed12

chrom - The name of the chromosome on which the genome feature exists.
    Any string can be used. For example, “chr1”, “III”, “myChrom”, “contig1112.23”.
    This column is required.
start - The zero-based starting position of the feature in the chromosome.
    The first base in a chromosome is numbered 0.
    The start position in each BED feature is therefore interpreted to be 1 greater than the start position listed in the feature. For example, start=9, end=20 is interpreted to span bases 10 through 20,inclusive.
    This column is required.
end - The one-based ending position of the feature in the chromosome.
    The end position in each BED feature is one-based. See example above.
    This column is required.
name - Defines the name of the BED feature.
    Any string can be used. For example, “LINE”, “Exon3”, “HWIEAS_0001:3:1:0:266#0/1”, or “my_Feature”.
    This column is optional.
score - The UCSC definition requires that a BED score range from 0 to 1000, inclusive. However, bedtools allows any string to be stored in this field in order to allow greater flexibility in annotation features. For example, strings allow scientific notation for p-values, mean enrichment values, etc. It should be noted that this flexibility could prevent such annotations from being correctly displayed on the UCSC browser.
    Any string can be used. For example, 7.31E-05 (p-value), 0.33456 (mean enrichment value), “up”, “down”, etc.
    This column is optional.
strand - Defines the strand - either ‘+’ or ‘-‘.
    This column is optional.
thickStart - The starting position at which the feature is drawn thickly.
    Allowed yet ignored by bedtools.
thickEnd - The ending position at which the feature is drawn thickly.
    Allowed yet ignored by bedtools.
itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0).
    Allowed yet ignored by bedtools.
blockCount - The number of blocks (exons) in the BED line.
    Allowed yet ignored by bedtools.
blockSizes - A comma-separated list of the block sizes.
blockStarts - A comma-separated list of block starts.
    Allowed yet ignored by bedtools.
"""
import argparse
import os
import re
import sys

from tqdm import tqdm

__author__ = "Zhang Yiming"
__since__ = "2020.01.13"


class Gtf2Bed12(object):

    def __init__(self):
        u"""
        init this class
        """
        args = self.argument_parser()
        print(args)
        self.input = os.path.abspath(args.input)
        self.output = os.path.abspath(args.output)
        self.check_dir()

        self.genes = {}
        self.transcripts = {}
        self.convert()
        pass

    def check_dir(self):
        u"""
        检查输入文件
        检查输出文件夹
        """
        if not os.path.isfile(self.input):
            raise FileNotFoundError("%s not found" % self.input)

        out_dir = os.path.dirname(self.output)

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    @staticmethod
    def argument_parser():
        u"""
        argument_parser
        """
        parser = argparse.ArgumentParser(
            description="Convert gtf|gff3 to bed12"
        )

        parser.add_argument(
            "-i",
            "--input",
            help="Path to input file",
            required=True
        )

        parser.add_argument(
            "-o",
            "--output",
            help="Path to output file",
            required=True
        )

        if len(sys.argv[1:]) <= 0:
            parser.print_help()
            exit(0)

        try:
            return parser.parse_args(sys.argv[1:])
        except argparse.ArgumentError as err:
            print(err, file=sys.stderr)
            parser.print_usage()
            exit(0)

    @staticmethod
    def __split_gtf_details__(line):
        u"""
        format gtf detailed messages into dict
        :param line: columns like gene_id "gene"; transcript_id "transcript";
        :return: dict {"gene_id": "gene"}
        """

        data = {}
        for message in line.strip().split(";"):
            if not message:
                continue

            key, value = re.split(r"[=\s]", message.strip())[:2]
            data[key.lower()] = re.sub(r"[\";]", "", value if ":" not in message else value.split(":")[1])
        return data

    @staticmethod
    def __get_value_from_data__(data, target, pop=True):
        u"""
        从gtf的详细列中构建出的字典，从其中提取出所需要的数据
        但是由于有多个可能性，比如：gene_id, geneID, ID等等，不同的标准下的gtf文件，太烦人了
        因此，在此通过正则来处理这个问题
        :param data: 从gtf文件中，提取出的字典信息
        :param target: 所要提取数据的目标，为gene_id等标签
        :return: string
        """
        ids = [target, target.replace("_", ""), target.split("_")[-1]]
        for i in ids:
            if i in data.keys():
                return data.pop(i) if pop else data[i]
        return "NA"

    @staticmethod
    def __format_bed12_exons__(data: (int, list)) -> str:
        u"""
        将外显子信息转换为三个值，总外显子数，外显子长度以及外显子的相对位置
        :param data: 由self.__split_gtf_details__构造的字典
        :return: string
        """
        data = sorted(data, key=lambda x: x[0])

        sizes = []
        pos = []
        for i in data:
            sizes.append(str(i[1] - i[0]))
            pos.append(str(i[0] - data[0][0]))

        return "\t".join([str(len(data)), ",".join(sizes), ",".join(pos)])

    def convert(self):
        u"""
        进行转化
        :return:
        """
        transcripts = {}
        with open(self.output, "w+") as w:
            with open(self.input) as r:
                for line in tqdm(r, desc="Reading"):
                    if line.startswith("#"):
                        continue

                    lines = line.split("\t")
                    data = self.__split_gtf_details__(lines[8])

                    if re.search("(transcript|mRNA)", lines[2], re.I):

                        transcript_id = self.__get_value_from_data__(data, "transcript_id", True)
                        if transcript_id == "NA":
                            transcript_id = self.__get_value_from_data__(data, "ID", True)

                        transcripts[transcript_id] = [
                            lines[0], lines[3], lines[4], transcript_id, "255", lines[6],
                            lines[3], lines[4], "255,0,0", []
                        ]

                    elif lines[2] == "exon":
                        parent = self.__get_value_from_data__(data, "transcript_id", False)
                        if parent == "NA":
                            parent = self.__get_value_from_data__(data, "Parent", False)

                        temp = transcripts.get(parent, [])
                        temp[-1].append([int(lines[3]), int(lines[4])])
                        transcripts[parent] = temp

            for transcript in tqdm(transcripts.values()):
                transcript[-1] = self.__format_bed12_exons__(transcript[-1])
                w.write("\t".join(transcript) + "\n")


if __name__ == '__main__':
    Gtf2Bed12()

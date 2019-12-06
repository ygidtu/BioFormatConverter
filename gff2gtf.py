#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
gff3 to gtf format
"""
import argparse
import os
import re
import sys

from tqdm import tqdm

__author__ = "Zhang Yiming"
__version__ = "0.1.1"
__since__ = 20180927


class Gff2Gtf(object):
    u"""
    convert gff3 to gtf
    """

    def __init__(self):
        u"""
        """
        args = self.argument_parser()
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

        outdir = os.path.dirname(self.output)

        if not os.path.exists(outdir):
            os.makedirs(outdir)

    @staticmethod
    def argument_parser():
        u"""
        argument_parser
        """
        parser = argparse.ArgumentParser(
            description="Convert gff3 to gtf"
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
    def split_gff_detail(line):
        u"""
        将gff中的信息列，拆分成字典
        """
        res = {}

        for i in line.split(";"):
            i = i.split("=")

            res[i[0]] = i[1] if ":" not in i[1] else i[1].split(":")[1]
        return res

    @staticmethod
    def concat_dict_to_string(data):
        u"""
        将字典中的键值对构成字符串
        """
        res = []

        for k, v in data.items():
            if "id" in k or "name" in k:
                res.append("%s \"%s\"" % (k, v))
        return "; ".join(res)

    def convert(self):
        u"""
        start to converting
        """
        with open(self.output, "w+") as w:
            w.write("#gtf-version")
            with open(self.input) as r:
                for line in tqdm(r):
                    if line.startswith("#"):
                        w.write(line)
                        continue

                    lines = line.rstrip().split("\t")
                    info = self.split_gff_detail(lines[-1])

                    # first class. eg: gene
                    if "ID" in info.keys() and \
                        "Parent" not in info.keys() and \
                            "gene" in lines[2]:

                        if "gene_id" not in info.keys():
                            info["gene_id"] = info["ID"]

                        if "gene_name" not in info.keys():
                            info["gene_name"] = info["Name"]

                        self.genes[info["gene_id"]] = info["gene_name"]

                    elif "Parent" in info.keys():

                        # second class. eg: transcripts
                        if info["Parent"] in self.genes.keys():
                            info["gene_id"] = info["Parent"]
                            info["gene_name"] = self.genes[
                                info["Parent"]
                            ]

                            if "transcript_id" not in info.keys():
                                info["transcript_id"] = info["ID"]

                            if "transcript_name" not in info.keys():
                                info["transcript_name"] = info["Name"] if "Name" in info.keys(
                                ) else "NA"

                            self.transcripts[info["transcript_id"]] = [
                                info["transcript_name"],
                                info["gene_id"]
                            ]

                            info["transcript_type"] = lines[2]
                            lines[2] = "transcript"

                        # third class. eg: exons
                        else:
                            if "transcript_id" not in info.keys():
                                info["transcript_id"] = info["Parent"]

                            transcript = self.transcripts[
                                info["Parent"]
                            ]

                            if "transcript_name" not in info.keys():
                                info["transcript_name"] = transcript[0]

                            if "gene_id" not in info.keys():
                                info["gene_id"] = transcript[1]

                            if "gene_name" not in info.keys():
                                info["gene_name"] = self.genes[
                                    transcript[1]
                                ]

                            ele_id = "%s_id" % lines[2].lower()
                            ele_name = "%s_name" % lines[2].lower()

                            if ele_id not in info.keys() and \
                                    "ID" in info.keys():
                                info[ele_id] = info.pop("ID")

                            if ele_name not in info.keys() and \
                                    "Name" in info.keys():
                                info[ele_name] = info.pop("Name")

                    # others just convert and output
                    else:
                        pass

                    lines[-1] = self.concat_dict_to_string(info)

                    if lines[2] in ("gene", "transcript", "exon"):
                        w.write("\t".join(lines) + "\n")


if __name__ == '__main__':
    Gff2Gtf()

#!/usr/bin/env python3
# -*- coding:utf-8 -*-

u"""
convert gtf files to gff3
"""
import argparse
import os
import re
import sys

__author__ = "Zhang Yiming"
__since__ = "2018.10.24"


class Gtf2Gff(object):

    def __init__(self):
        u"""
        init this class
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

        out_dir = os.path.dirname(self.output)

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    @staticmethod
    def argument_parser():
        u"""
        argument_parser
        """
        parser = argparse.ArgumentParser(
            description="Convert gtf to gff3"
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
            key, value = message.strip().split(" ")
            data[key.lower()] = re.sub(r"[\";]", "", value)
        return data

    @staticmethod
    def __get_value_from_data__(data, target):
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
                return data.pop(i)
        return "NA"

    def __format_gff_details__(self, data, label):
        u"""
        将获取到的gtf的信息，format成gff3样式
        :param data: 由self.__split_gtf_details__构造的字典
        :param label: gtf文件，第二列表明的元件类型
        :return: string
        """
        result = ""
        if label == "gene":
            result += "ID=%s;Name=%s" % (
                self.__get_value_from_data__(data=data, target="gene_id"),
                self.__get_value_from_data__(data=data, target="gene_name"),
            )
        elif label == "transcript":
            result += "ID=%s;Name=%s;Parent=%s" % (
                self.__get_value_from_data__(data=data, target="transcript_id"),
                self.__get_value_from_data__(data=data, target="transcript_name"),
                self.__get_value_from_data__(data=data, target="gene_id"),
            )

            if "gene_name" in data.keys():
                data.pop("gene_name")
        else:
            ids = ["exon_id", "protein_id", "%s_id" % label.lower(), "ID"]

            for i in ids:
                if i in data.keys():
                    ids = data.pop(i)
                    break

            ids = ids if isinstance(ids, str) else "NA"

            parent = self.__get_value_from_data__(data=data, target="transcript_id")

            if "exon_number" in data.keys() and ids == "NA":
                ids = "%s.%s" % (parent, data["exon_number"])

            result += "ID=%s;Parent=%s" % (ids, parent)

            if "gene_name" in data.keys():
                data.pop("gene_name")

        for key, value in data.items():
            result += ";%s=%s" % (key, value)

        return result

    def convert(self):
        u"""
        进行转化
        :return:
        """
        with open(self.output, "w+") as w:
            with open(self.input) as r:
                for line in r:
                    if line.startswith("#"):
                        continue

                    lines = line.split("\t")
                    data = self.__split_gtf_details__(lines[8])
                    lines[8] = self.__format_gff_details__(data, lines[2])

                    w.write("\t".join(lines) + "\n")


if __name__ == '__main__':
    Gtf2Gff()

#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
从gmap的splicesite中提取位点
"""
import re

from fire import Fire
from tqdm import tqdm


class converter(object):

    def __init__(self, infile, outfile):
        self.infile = infile
        self.outfile = outfile

        self.data = []
        self.junctions = {}
        self.__convert__()

    def __convert__(self):
        current = []
        strand = "."
        chromosome = "."
        pattern = r"^\s+(?P<strand>[+-])(?P<chrom>[\w\.]+):(?P<start>\d+)-(?P<end>\d+)\s+\(\d+-\d+\)\s+\d+%.*"
        with open(self.infile) as r:
            for line in tqdm(r):
                data = re.search(pattern, line)

                if not data:
                    if current:
                        current = sorted(current)

                        gene_range = "%s\t%d\t%d\t%s" % (
                            chromosome,
                            current[0],
                            current[-1],
                            strand
                        )

                        junctions = []
                        for i in range(1, len(current) - 1, 2):
                            junctions.append(current[i] + 1)
                            junctions.append(current[i + 1] - 1)

                            key = "%s\t%d\t%d\t%s" % (
                                chromosome,
                                current[i] + 1,
                                current[i + 1] - 1,
                                strand
                            )

                            if key in self.junctions.keys():
                                self.junctions[key] += 1
                            else:
                                self.junctions[key] = 1

                        current.clear()

                        self.data.append("%s\t%s" % (
                            gene_range, ",".join([str(x) for x in junctions])
                        ))
                    continue

                chromosome, strand = data["chrom"], data["strand"]
                current.append(int(data["start"]))
                current.append(int(data["end"]))

        with open(self.outfile, "w+") as w:
            w.write("\n".join(self.data))

        with open(self.outfile + ".junctions", "w+") as w:
            for key, value in self.junctions.items():
                w.write("%s\t%d\n" % (key, value))


if __name__ == '__main__':
    Fire(converter)

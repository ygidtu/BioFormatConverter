# BioFormatConverter
一些常见生信文件互相转化工具，发现需求还挺大的，但是现有工具各种有问题，就自己现需要先写吧

- gff2gtf.py -> gff3格式转化为将诶gtf格式。已在ensembl和NCBI的格式上进行过测试
- gmap_splicesites2sj.py -> gmap -A输出的alignment情况，提取出两个文件，一个包含reads位点和intron sites；另一个包含junctions的位点和count
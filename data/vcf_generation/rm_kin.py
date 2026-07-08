#! /usr/bin/env python3

with open("kin_to_rm.txt") as f:
    kinrm = f.read().splitlines()

sentinel = False

with open("YFT.pcrelate.kinrm.haplo.gen") as f:
    with open("CORRECTYFT.pcrelate.kinrm.haplo.gen", "w") as o:
        for i in f:
            for k in kinrm:
                if i.startswith(k):
                    sentinel = True
                    break
            if not sentinel:
                o.write(i)
            sentinel = False

sentinel = False

with open("YFT.snp.kinrm.pcrelate.gen") as f:
    with open("CORRECTYFT.YFT.snp.kinrm.pcrelate.gen", "w") as o:
        for i in f:
            for k in kinrm:
                if i.startswith(k):
                    sentinel = True
                    break
            if not sentinel:
                o.write(i)
            sentinel = False

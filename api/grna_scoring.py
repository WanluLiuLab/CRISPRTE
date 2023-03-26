# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#
# gRNA scoring for 20bp gRNA sequence
# 1. Moreno-Mateos, M. A. et al. CRISPRscan: designing highly efficient sgRNAs for CRISPR-Cas9 targeting in vivo. Nat Methods 12, 982–988 (2015).
# 2. John G. Doench*, Nicolo Fusi*, Meagan Sullender*, Mudra Hegde*, Emma W. Vaimberg*, Katherine F. Donovan, Ian Smith, Zuzana Tothova, Craig Wilen , Robert Orchard, Herbert W. Virgin, Jennifer Listgarten*, David E. Root. Optimized sgRNA design to maximize activity and minimize off-target effects for genetic screens with CRISPR-Cas9. Nature Biotechnology, 2016. 

# ------------------------- #
# Python Modules
# ------------------------- #

import numpy as np
from azimuth.model_comparison import predict
from typing import List

def calculate_azimuth(sequences:list):
	return predict(np.array(sequences))




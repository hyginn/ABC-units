# 2021-11-30_In-Class-SARS-COV-2-Omicron
#
#         SARS COV-2 Omicron variant
#         Annotation example
#
# Code and comments for BCH441 in-class exploration, Tuesday, 2021-11-30
# Explorers:  N/A
# Scribe:     boris.steipe@utoronto.ca
#
# ==============================================================================
#
#
# 1. How does does SARS-COV-2 work?
# 1.1  capsid and spike protein

# Fetch an article from https://pubmed.ncbi.nlm.nih.gov/32694201/
# Cai et al.
# "Distinct conformational states of SARS-CoV-2 spike protein"
# Science. 2020 Sep 25;369(6511):1586-1592
#   Image:
#

# Load things in Chimera ...
CXPORT <- 61803

CX("open 6XR8")
CX("camera sbs")
CX("color /A #8F89F7")
CX("color /B #996F6F")
CX("color /C #AF5200")
CX("hide :NAG,MAN,FUC")
CX("color /B,C #32384F")


# === receptor binding fragment (S1)  =======

CX("select /A:14-305") # "NTD"
CX("cofr sel")
CX("color sel #89a3f7")

CX("select /A:306-330") # connecting peptide
CX("cofr sel")
CX("color sel #999999")

CX("select /A:331-528") # RBD
CX("cofr sel")
CX("color sel #89d2f7")

CX("select /A:323-330,529-591") # CTD1
CX("cofr sel")
CX("color sel #89f7cb")

CX("select /A:676-689") # Furin cleavage site in the disordered loop
                        # sequence (T)qtnspRRARsva(S)
CX("cofr sel")
CX("color sel #f789a6")

CX("select /A:309-315,592-698") # CTD2 (technically, after 686 - furin
                                # cleavage site it is considered S2)
CX("cofr sel")
CX("color sel #b1f789")

# === fusion fragment (S2)  =======
CX("select /A:699-910") # S2
CX("cofr sel")
CX("color sel #e8f789")


CX("select /A:816") # second cleavage site (TMPRSS2 target)
CX("cofr sel")
CX("color sel #9df789")

CX("select /A:817-834") # FP FUSION PEPTIDE
CX("cofr sel")
CX("color sel #f3f789")

CX("select /A:817-834") # FPPR fusion peptide proximal region
CX("cofr sel")
CX("color sel #f7e589")

CX("select /A:910-985") # HR1 Heptad repeat region 1
CX("cofr sel")
CX("color sel #f7cb89")

CX("select /A:986-1035") # CH coiled-coil core
CX("cofr sel")
CX("color sel #f7b989")

CX("select /A:1035-1068") # CD connector domain
CX("cofr sel")
CX("color sel #f78589")

CX("select /A:1063-1211") # HR2 Heptad repeat 2
CX("cofr sel")
CX("color sel #f78d89")

# Omicron mutations 6M0J

CX("select /A:67,95,142,212,339,371,373,375,417,440,446,477,478,484,493,496,498,501,505,547,614,655,679,681,764,796,856,954,969,981")
CX("color sel #fd6328")
CX("cofr /A:67")
CX("cofr /A:440")

CX("color /A #32384F")

CX("select /B:67,95,142,212,339,371,373,375,417,440,446,477,478,484,493,496,498,501,505,547,614,655,679,681,764,796,856,954,969,981")
CX("color sel #fd6328")

# ACE2 complex 6M0J
#


# [END]

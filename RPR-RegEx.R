# RPR-RegEx.R
#
# Purpose: A Bioinformatics Course:
#              R code accompanying the RPR-RegEx unit
#
# Version: 0.1
#
# Date:    2017  08  25
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 0.1    First code
#
# TODO:
#
#
# == HOW TO WORK WITH LEARNING UNIT FILES ======================================
#
# DO NOT SIMPLY  source()  THESE FILES!

# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
#  going on. That's not how it works ...
#
# ==============================================================================

# ==============================================================================
#        PART ONE: REGEX EXAMPLE
# ==============================================================================

# The Mbp1 sequence as copied from the NCBI Website
mySeq <- "
1 msnqiysary sgvdvyefih stgsimkrkk ddwvnathil kaanfakakr trilekevlk
61 ethekvqggf gkyqgtwvpl niakqlaekf svydqlkplf dftqtdgsas pppapkhhha
121 skvdrkkair sastsaimet krnnkkaeen qfqsskilgn ptaaprkrgr pvgstrgsrr
181 klgvnlqrsq sdmgfprpai pnssisttql psirstmgpq sptlgileee rhdsrqqqpq
241 qnnsaqfkei dledglssdv epsqqlqqvf nqntgfvpqq qssliqtqqt esmatsvsss
301 pslptspgdf adsnpfeerf pgggtspiis miprypvtsr pqtsdindkv nkylsklvdy
361 fisnemksnk slpqvllhpp phsapyidap idpelhtafh wacsmgnlpi aealyeagts
421 irstnsqgqt plmrsslfhn sytrrtfpri fqllhetvfd idsqsqtvih hivkrksttp
481 savyyldvvl skikdfspqy rielllntqd kngdtalhia skngdvvffn tlvkmgaltt
541 isnkegltan eimnqqyeqm miqngtnqhv nssntdlnih vntnnietkn dvnsmvimsp
601 vspsdyityp sqiatnisrn ipnvvnsmkq masiyndlhe qhdneikslq ktlksisktk
661 iqvslktlev lkesskdeng eaqtnddfei lsrlqeqntk klrkrliryk rlikqkleyr
721 qtvllnklie detqattnnt vekdnntler lelaqeltml qlqrknklss lvkkfednak
781 ihkyrriire gtemnieevd ssldvilqtl iannnknkga eqiitisnan sha
//
"

mySeq     # "\n" means: line-break

mySeq <- gsub("[^a-zA-Z]", "", mySeq) # replace all non-letters with ""

mySeq

# Now to change the sequence to upper-case. R has toupper()
# and tolower().

toupper(mySeq)

# CONTINUE ...






# [END]

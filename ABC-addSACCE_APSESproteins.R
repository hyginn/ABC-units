# addSACCE_APSESproteins.R
# Adds the Saccharomyces cerevisiae APSES proteins to myDB
#

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbAutoincrement(myDB$protein$ID, ns = "ref"),
              name = "SWI4_SACCE",
              RefSeqID = "NP_011036",
              UniProtID = "P25302",
              taxonomy.ID = as.integer(4932),
              sequence = dbSanitizeSequence("
        1 mpfdvlisnq kdntnhqnit pisksvllap hsnhpvieia tysetdvyec yirgfetkiv
       61 mrrtkddwin itqvfkiaqf sktkrtkile kesndmqhek vqggygrfqg twipldsakf
       121 lvnkyeiidp vvnsiltfqf dpnnpppkrs knsilrktsp gtkitspssy nktprkknss
       181 sstsatttaa nkkgkknasi nqpnpsplqn lvfqtpqqfq vnssmnimnn ndnhttmnfn
       241 ndtrhnlinn isnnsnqsti iqqqksihen sfnnnysatq kplqffpipt nlqnknvaln
       301 npnnndsnsy shnidnvins snnnnngnnn nliivpdgpm qsqqqqqhhh eyltnnfnhs
       361 mmdsitngns kkrrkklnqs neqqfynqqe kiqrhfklmk qpllwqsfqn pndhhneycd
       421 sngsnnnnnt vasngssiev fssnendnsm nmssrsmtpf sagntssqnk lenkmtdqey
       481 kqtiltilss erssdvdqal latlypapkn fninfeiddq ghtplhwata maniplikml
       541 itlnanalqc nklgfncitk sifynncyke nafdeiisil kiclitpdvn grlpfhylie
       601 lsvnksknpm iiksymdsii lslgqqdynl lkiclnyqdn igntplhlsa lnlnfevynr
       661 lvylgastdi lnldnespas imnkfntpag gsnsrnnntk adrklarnlp qknyyqqqqq
       721 qqqpqnnvki pkiiktqhpd kedstadvni aktdsevnes qylhsnqpns tnmntimedl
       781 sninsfvtss vikdikstps kilenspily rrrsqsisde kekakdnenq vekkkdplns
       841 vktampsles pssllpiqms plgkyskpls qqinklntkv sslqrimgee iknldnevve
       901 tessisnnkk rlitiahqie dafdsvsnkt pinsisdlqs riketsskln sekqnfiqsl
       961 eksqalklat ivqdeeskvd mntnssshpe kqedeepipk stsetsspkn tkadakfsnt
       1021 vqesydvnet lrlateltil qfkrrmttlk iseakskins svkldkyrnl igitienids
       1081 klddiekdlr ana"),
              stringsAsFactors = FALSE))

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbAutoincrement(myDB$protein$ID, ns = "ref"),
              name = "PHD1_SACCE",
              RefSeqID = "NP_012881",
              UniProtID = "P36093",
              taxonomy.ID = as.integer(4932),
              sequence = dbSanitizeSequence("
        1 myhvpemrlh yplvntqsna aitptrsydn tlpsfnelsh qstinlpfvq retpnayanv
       61 aqlatsptqa ksgyycryya vpfptypqqp qspyqqavlp yatipnsnfq pssfpvmavm
      121 ppevqfdgsf lntlhphtel ppiiqntndt svarpnnlks iaaasptvta ttrtpgvsst
      181 svlkprvitt mwedenticy qveangisvv rradnnming tkllnvtkmt rgrrdgilrs
      241 ekvrevvkig smhlkgvwip ferayilaqr eqildhlypl fvkdiesivd arkpsnkasl
      301 tpksspapik qepsdnkhei ateikpksid alsngastqg agelphlkin hidteaqtsr
      361 aknels"),
              stringsAsFactors = FALSE))

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbAutoincrement(myDB$protein$ID, ns = "ref"),
              name = "SOK2_SACCE",
              RefSeqID = "NP_013729",
              UniProtID = "P53438",
              taxonomy.ID = as.integer(4932),
              sequence = dbSanitizeSequence("
        1 mpignpintn diksnrmrqe snmsavsnse stigqstqqq qqqqqylgqs vqplmpvsyq
       61 yvvpeqwpyp qyyqqpqsqs qqqlqsqpqm yqvqesfqss gsdsnasnpp stsvgvpsna
      121 tatalpngsa ittkksnnst nisnnvpyyy yfpqmqaqqs maysypqayy yypangdgtt
      181 ngatpsvtsn qvqnpnlekt ystfeqqqqh qqqqqlqaqt ypaqppkign afskfsksgp
      241 psdsssgsms pnsnrtsrns nsisslaqqp pmsnypqpst yqypgfhkts sipnshspip
      301 prslttptqg ptsqngplsy nlpqvgllpp qqqqqvsply dgnsitppvk pstdqetylt
      361 anrhgvsdqq ydsmaktmns fqtttirhpm pliattnatg sntsgtsasi irprvtttmw
      421 edektlcyqv eangisvvrr adndmvngtk llnvtkmtrg rrdgilkaek irhvvkigsm
      481 hlkgvwipfe ralaiaqrek iadylyplfi rdiqsvlkqn npsndsssss sstgiksisp
      541 rtyyqpinny qnpngpsnis aaqltyssmn lnnkiipnns ipavstiaag ekplkkctmp
      601 nsnqleghti tnlqtlsatm pmkqqlmgni asplsyprna tmnsastlgi tpadskpltp
      661 sptttntnqs sesnvgsiht gitlprvese sashskwske adsgntvpdn qtlkeprssq
      721 lpisaltstd tdkiktstsd eatqpnepse aepvkesess ksqvdgagdv sneeiaaddt
      781 kkqek"),
              stringsAsFactors = FALSE))

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbAutoincrement(myDB$protein$ID, ns = "ref"),
              name = "XBP1_SACCE",
              RefSeqID = "NP_012165",
              UniProtID = "P40489",
              taxonomy.ID = as.integer(4932),
              sequence = dbSanitizeSequence("
        1 mkypafsins dtvhltdnpl ddyqrlylvs vldrdsppas fsaglnirkv nykssiaaqf
       61 thpnfiisar dagngeeaaa qnvlncfeyq fpnlqtiqsl vheqtllsql assatphsal
      121 hlhdknilmg kiilpsrsnk tpvsasptkq ekkalstasr enatssltkn qqfkltkmdh
      181 nlindklinp nncviwshds gyvfmtgiwr lyqdvmkgli nlprgdsvst sqqqffckae
      241 fekilsfcfy nhssftsees ssvllsssts sppkrrtstg stfldanass sstsstqann
      301 yidfhwnnik pelrdlicqs ykdflinelg pdqidlpnln panftkrirg gyikiqgtwl
      361 pmeisrllcl rfcfpiryfl vpifgpdfpk dceswylahq nvtfassttg agaataataa
      421 antstnftst avarprqkpr prprqrstsm shskaqklvi edalpsfdsf venlglssnd
      481 knfikknskr qksstytsqt sspigprdpt vqilsnlasf ynthghrysy pgniyipqqr
      541 yslpppnqls spqrqlnyty dhihpvpsqy qsprhynvps spiapapptf pqpygddhyh
      601 flkyasevyk qqnqrpahnt ntnmdtsfsp rannslnnfk fktnskq"),
              stringsAsFactors = FALSE))

# [END]

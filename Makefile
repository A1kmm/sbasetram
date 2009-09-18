all:
	ghc -O2 -Wall -fno-warn-missing-signatures -fvia-c -optc-O3 --make gmatim.hs
	ghc -O2 -Wall -fno-warn-missing-signatures -fvia-c -optc-O3 --make sbasetram.hs
	ghc -O2 -Wall -fno-warn-missing-signatures -fvia-c -optc-O3 --make pfilter.hs
	ghc -O2 -Wall -fno-warn-missing-signatures -fvia-c -optc-O3 --make BuildGeneList.hs

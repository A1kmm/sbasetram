all:
	ghc -O2 -Wall -fno-warn-missing-signatures -fvia-c -optc-O3 --make gmatim.hs

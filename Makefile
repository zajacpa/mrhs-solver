mrhs.hillc: mrhs.bm.c mrhs.bv.c mrhs.c mrhs.hillc.c mrhs.rz.c mrhs.tester.c mrhs.1.7.c
	gcc -o mrhs mrhs.bm.c mrhs.bv.c mrhs.c mrhs.hillc.c mrhs.rz.c mrhs.tester.c mrhs.1.7.c -I . -lm -D_VERBOSITY=4

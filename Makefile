mrhs.hillc: mrhs.bm.c mrhs.c mrhs.hillc.0.2.c mrhs.hillc.tester.c
	gcc -o mrhs mrhs.bm.c mrhs.c mrhs.hillc.0.2.c mrhs.hillc.tester.c -I . -lm -D_VERBOSITY=4

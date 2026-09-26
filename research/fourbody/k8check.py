exec(open("k7check.py").read().split("import sys\nKK")[0])
for K in [16, 20]:
    print("K=%d: rhombic %.5f   hex %.5f   square %.5f" % (K, Qk(0.5, 1.306706, 8, K), Qk(0.5, 0.8660254, 8, K), Qk(0, 1, 8, K)), flush=True)

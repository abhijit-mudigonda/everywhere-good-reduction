if __name__ == "__main__":
    with open("ub-annoying_d.txt", 'r') as goodfile, open("tmp-annoying.txt",'w') as fout:
        for line in goodfile:
            d = Integer(line)
            if d.is_squarefree():
                fout.write(str(d)+'\n')

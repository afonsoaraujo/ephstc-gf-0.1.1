import numpy, sys, os, time, glob
import matplotlib.pyplot as plt
##############################################################################################
ncol=63;nlin=36 # dimensao da matriz
os.chdir('/home/afonsoaraujo/work/bdshared/negro/forc/gldas/lwdown/') # define diretorio de arquivos de saida
################################################################################################
t = [] # criando uma matriz aonde vou inserir meus dados
file = sorted(glob.glob('/home/afonsoaraujo/work/bdshared/negro/forc/gldas/lwdown/*.r4')) # funcao que lista todos os arquivos (*) da pasta referida e retorna um vetor de nomes  
#nk = len(file)
nk = 2920
c = [ [ 0 for nlnt in range(ncol*nlin) ] for nt in range(nk) ]

fo = open('stations.txt', r+)

for k in range (0,2):
    a = file[k]
    tm = len(file[k])
    b =  numpy.fromfile(a, dtype=numpy.float32)
    for i in range (0,ncol):
       for j in range (0,nlin):
          l = i*j
#          print k, i, j, l, b[l]
          c[l][k] = b[l]
    print c[l][k], l, k
#   fo.write(c[l][k]\n, fmt "%s")
  
#numpy.savetxt('stations.txt', c[l][k], delimiter = ",", fmt = "%s")


#for i in range (0,ncol):
#   for j in range (0,nlin):
#      l = i*j
#      print k,i,j,l, c[k][l]
#      plt.plot(k,c[k][l])
#      plt.axis([0, nk, -100, 500])        
#      plt.show()	  

       



##################################################################################################

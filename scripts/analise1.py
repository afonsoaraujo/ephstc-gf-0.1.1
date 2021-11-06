import numpy, sys, os, time, glob
import matplotlib.pyplot as plt
##############################################################################################
ncol=63;nlin=36 # dimensao da matriz
os.chdir('/home/afonsoaraujo/work/bdshared/negro/obs/gldas/xteste/') # define diretorio de arquivos de saida
################################################################################################
t = [] # criando uma matriz aonde vou inserir meus dados
file = sorted(glob.glob('/home/afonsoaraujo/work/bdshared/negro/obs/gldas/xteste/*.r4')) # funcao que lista todos os arquivos (*) da pasta referida e retorna um vetor de nomes  
#ntemp = len(file)
ntemp = 2
c = [ [ [ 0 for nk in range(ntemp) ] for ni in range(ncol) ] for nj in range(nlin) ]

for j in range (0,nlin):
   for i in range (0,ncol):
      for k in range (0,ntemp): 
         a  = file[k]
         b  = numpy.fromfile(a, dtype=numpy.float32)
	 l = j + i
	 c[j][i][k] = b[l]
#	 print k, c[j][i][k]


for j in range (0,nlin):
   for i in range (0,ncol):
      print k, c[:][:][k] 
   
#       print c[:][i][j]
#      ioff()
#      plt.plot(c[:][i][j])
#      plt.show()	  

       



##################################################################################################

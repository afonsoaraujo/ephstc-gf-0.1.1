import numpy, sys, os, time, glob
##############################################################################################
m=63;n=36 # dimensao da matriz
os.chdir('C:\\Users\\Jmarcos\\Desktop\\negro\\forc\\gldas') # define diretorio de arquivos de saida
################################################################################################
t = [] # criando uma matriz aonde vou inserir meus dados
file = sorted(glob.glob('C:\\Users\\Jmarcos\\Desktop\\negro\\forc\\gldas\\lwdown\\*')) # funcao que lista todos os arquivos (*) da pasta referida e retorna um vetor de nomes  
for i in range(len(file)): # len file define o tamanho na matriz e faco loop comecando com o elemento i=0 ateh o ultimo elemento da matriz
    a = file[i]
    tm = len(file[i])
    b =  numpy.fromfile(a, dtype=numpy.float32)
    mask = numpy.not_equal(b,9.999e+20)
    b = mask*b
    media = numpy.average(b);dpad  = numpy.std(b);max   = numpy.amax(b);min   = numpy.amin(b)
    t = numpy.append(t,(a[(tm-11):(tm-3)], media,dpad,max,min))
t = numpy.reshape(t,(len(t)/5,5))
numpy.savetxt('lwdown.txt', t, delimiter = ",", fmt = "%s")
##################################################################################################
t = [] # criando uma matriz aonde vou inserir meus dados
file = sorted(glob.glob('C:\\Users\\Jmarcos\\Desktop\\negro\\forc\\gldas\\psurf\\*')) # funcao que lista todos os arquivos (*) da pasta referida e retorna um vetor de nomes  
for i in range(len(file)): # len file define o tamanho na matriz e faco loop comecando com o elemento i=0 ateh o ultimo elemento da matriz
    a = file[i]
    tm = len(file[i])
    b =  numpy.fromfile(a, dtype=numpy.float32)
    mask = numpy.not_equal(b,9.999e+20)
    b = mask*b
    media = numpy.average(b);dpad  = numpy.std(b);max   = numpy.amax(b);min   = numpy.amin(b)
    t = numpy.append(t,(a[(tm-11):(tm-3)], media,dpad,max,min))
t = numpy.reshape(t,(len(t)/5,5))
numpy.savetxt('psurf.txt', t, delimiter = ",", fmt = "%s")
##############################################################################################
t = [] # criando uma matriz aonde vou inserir meus dados
file = sorted(glob.glob('C:\\Users\\Jmarcos\\Desktop\\negro\\forc\\gldas\\qair\\*')) # funcao que lista todos os arquivos (*) da pasta referida e retorna um vetor de nomes  
for i in range(len(file)): # len file define o tamanho na matriz e faco loop comecando com o elemento i=0 ateh o ultimo elemento da matriz
    a = file[i]
    tm = len(file[i])
    b =  numpy.fromfile(a, dtype=numpy.float32)
    mask = numpy.not_equal(b,9.999e+20)
    b = mask*b
    media = numpy.average(b);dpad  = numpy.std(b);max   = numpy.amax(b);min   = numpy.amin(b)
    t = numpy.append(t,(a[(tm-11):(tm-3)], media,dpad,max,min))
t = numpy.reshape(t,(len(t)/5,5))
numpy.savetxt('qair.txt', t, delimiter = ",", fmt = "%s")
##############################################################################################
t = [] # criando uma matriz aonde vou inserir meus dados
file = sorted(glob.glob('C:\\Users\\Jmarcos\\Desktop\\negro\\forc\\gldas\\rainf\\*')) # funcao que lista todos os arquivos (*) da pasta referida e retorna um vetor de nomes  
for i in range(len(file)): # len file define o tamanho na matriz e faco loop comecando com o elemento i=0 ateh o ultimo elemento da matriz
    a = file[i]
    tm = len(file[i])
    b =  numpy.fromfile(a, dtype=numpy.float32)
    mask = numpy.not_equal(b,9.999e+20)
    b = mask*b
    media = numpy.average(b);dpad  = numpy.std(b);max   = numpy.amax(b);min   = numpy.amin(b)
    t = numpy.append(t,(a[(tm-11):(tm-3)], media,dpad,max,min))
t = numpy.reshape(t,(len(t)/5,5))
numpy.savetxt('rainf.txt', t, delimiter = ",", fmt = "%s")
################################################################################################
t = [] # criando uma matriz aonde vou inserir meus dados
file = sorted(glob.glob('C:\\Users\\Jmarcos\\Desktop\\negro\\forc\\gldas\\swdown\\*')) # funcao que lista todos os arquivos (*) da pasta referida e retorna um vetor de nomes  
for i in range(len(file)): # len file define o tamanho na matriz e faco loop comecando com o elemento i=0 ateh o ultimo elemento da matriz
    a = file[i]
    tm = len(file[i])
    b =  numpy.fromfile(a, dtype=numpy.float32)
    mask = numpy.not_equal(b,9.999e+20)
    b = mask*b
    media = numpy.average(b);dpad  = numpy.std(b);max   = numpy.amax(b);min   = numpy.amin(b)
    t = numpy.append(t,(a[(tm-11):(tm-3)], media,dpad,max,min))
t = numpy.reshape(t,(len(t)/5,5))
numpy.savetxt('swdown.txt', t, delimiter = ",", fmt = "%s")
################################################################################################
t = [] # criando uma matriz aonde vou inserir meus dados
file = sorted(glob.glob('C:\\Users\\Jmarcos\\Desktop\\negro\\forc\\gldas\\tair\\*')) # funcao que lista todos os arquivos (*) da pasta referida e retorna um vetor de nomes  
for i in range(len(file)): # len file define o tamanho na matriz e faco loop comecando com o elemento i=0 ateh o ultimo elemento da matriz
    a = file[i]
    tm = len(file[i])
    b =  numpy.fromfile(a, dtype=numpy.float32)
    mask = numpy.not_equal(b,9.999e+20)
    b = mask*b
    media = numpy.average(b);dpad  = numpy.std(b);max   = numpy.amax(b);min   = numpy.amin(b)
    t = numpy.append(t,(a[(tm-11):(tm-3)], media,dpad,max,min))
t = numpy.reshape(t,(len(t)/5,5))
numpy.savetxt('tair.txt', t, delimiter = ",", fmt = "%s")
################################################################################################
t = [] # criando uma matriz aonde vou inserir meus dados
file = sorted(glob.glob('C:\\Users\\Jmarcos\\Desktop\\negro\\forc\\gldas\\wind\\*')) # funcao que lista todos os arquivos (*) da pasta referida e retorna um vetor de nomes  
for i in range(len(file)): # len file define o tamanho na matriz e faco loop comecando com o elemento i=0 ateh o ultimo elemento da matriz
    a = file[i]
    tm = len(file[i])
    b =  numpy.fromfile(a, dtype=numpy.float32)
    mask = numpy.not_equal(b,9.999e+20)
    b = mask*b
    media = numpy.average(b);dpad  = numpy.std(b);max   = numpy.amax(b);min   = numpy.amin(b)
    t = numpy.append(t,(a[(tm-11):(tm-3)], media,dpad,max,min))
t = numpy.reshape(t,(len(t)/5,5))
numpy.savetxt('wind.txt', t, delimiter = ",", fmt = "%s")
